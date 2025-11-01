#-------------------------------------------------------------------------------
# ==============================================================================
# This R script demonstrates how to apply PULSE method to correct LE...
# measurements for energy imbalance at eddy covariance flux sites.
# Author: Pushpendra Raghav (ppushpendra@ua.edu; University of Alabama)
# Oct 30, 2025
# ==============================================================================
# Some important functions
# ==============================================================================
is_outlier <- function(y, multiplier = 2) {
  abs(y - median(y, na.rm = TRUE)) > multiplier * sd(y, na.rm = TRUE)
}

quantile_reg <- function(x, y , PolyDeg=1, tau=0.95, weights){
  "
  Quantile regression
    Fits a polynomial function (of degree PolyDeg) using quantile regression based on a percentile (tau).
    Based on script by Dr. Phillip M. Feldman, and based on method by Koenker, Roger, and
    Gilbert Bassett Jr. Regression Quantiles. Econometrica: Journal of the Econometric Society, 1978, 33-50.
  
  Parameters
    ----------
    x : independent variable
    y : dependent variable
    PolyDeg : Degree of polynomial function
    tau : Percentile for the data to fit to [0-1]
    weights : Vector to weight each point, must be same size as x
    
  Returns
    -------
    The resulting parameters in order of degree from high to low
  
  "
  model <- function(x, beta){
    "
    This example defines the model as a polynomial, where the coefficients of the
    polynomial are passed via `beta`.
    "
    library(signal)
    if(PolyDeg==0){
      return(x*beta)
    } else{
      return(polyval(beta, x))
    }
  }
  N_coefficients <- PolyDeg+1
  
  tilted_abs <- function(tau, x, weights){
    "
     The tilted absolute value function is used in quantile regression.
      INPUTS
       tau: This parameter is a probability, and thus takes values between 0 and 1.
       x: This parameter represents a value of the independent variable, and in
       general takes any real value (float).
    "
    return (weights * x * (tau - (x < 0)))
  }
  objective <- function(beta, tau, weights){
    "
    The objective function to be minimized is the sum of the tilted absolute
    values of the differences between the observations and the model.
    "
    return(sum(tilted_abs(tau, y - model(x, beta), weights)))
  }
  # Build weights if they don't exits:
  if(is.na(weights)){
    weights <- rep(1,length(x))
  }
  # Define starting point for optimization:
  beta_0 <- rep(0, N_coefficients)  
  if(N_coefficients >= 2){
    beta_0[1] <- 1.0
  }
  # `beta_hat[i]` will store the parameter estimates for the quantile
  # corresponding to `fractions[i]`:
  library(stats)
  beta_hat <- optim(beta_0, objective,weights=weights, tau=tau)
  beta_hat <- beta_hat$par
  return(beta_hat)
}
# ==============================================================================
# Load or Install Required R packages
load_or_install <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
packages <- c(
  "lubridate",
  "dplyr",
  "data.table",
  "quantreg",
  "bigleaf"
)
suppressMessages(sapply(packages, load_or_install))
# ==============================================================================
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Data Needed: GPP, VPD, LE, Ta, NETRAD, G, H, and other supporting variables (if available)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Read Data from standard flux data (FLUXNET2015) file 
# (Note: I have extracted some key variables only to reduce file size. 
# However, this script will work with any FLUXNET2015 (or related) files)

# Set working directory
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  current_path <- rstudioapi::getActiveDocumentContext()$path
  setwd(dirname(current_path))
  cat("Working directory set to script folder:\n", getwd(), "\n")
} else {
  cat("Not running in RStudio. Working directory not changed.\n")
}

# Main
ds <- read.csv("sample_data_FR-Pue_2000_2014.csv")
ds[ds == -9999] <- NA
ds$DateTime <- as.POSIXct(as.character(ds$TIMESTAMP_START), format = "%Y%m%d%H%M", tz = "UTC")
ds$Year <- year(ds$DateTime)
# ----Data screening and quality control for the reference uWUEp estimation----
# --> (1) Use only half hourly data with high confidence, i.e., original or most reliable data acc. to quality flags (Quality mask)
# --> (2) only daylight data (i.e., when SW_IN_POT > 0) 
# --> (3) Select data during the growing season: i.e., data for days when when average half-hourly 
# GPP was at least 10% of the 95th percentile of all the half-hourly GPP for the site (SeasonMask)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# 1. ----Quality Mask of GPP data----
ds$QualityMask <- TRUE
ds$QualityMask[ds$NEE_VUT_USTAR50_QC > 2 | ds$NEE_VUT_USTAR50_QC < 0 | !is.finite(ds$GPP_NT_VUT_REF)] <- FALSE

# 2. ----Zero Mask ----
ds$ZeroMask <- TRUE
ds$ZeroMask[ds$SW_IN_POT <= 0] <- FALSE

# 3. ----Season Mask----
ds$SeasonMask <- TRUE
ds <- ds %>% 
  group_by(Date = date(DateTime)) %>% 
  mutate(GPP_day = mean(GPP_NT_VUT_REF, na.rm=T))
ds$SeasonMask[ds$GPP_day <= 0.10*quantile(ds$GPP_day,0.95, na.rm=T)] <- FALSE

ds <- ds %>% as.data.table()        # to simplify operations on ds    
#-------------------------------------------------------------------------------
ds$uWUEp_Mask   = ds$ZeroMask & ds$QualityMask & ds$SeasonMask
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------
ds$GPP <- ds$GPP_NT_VUT_REF*12.001/1e6 * 24*3600 # GPP from umolC m-2 s-1 to gC m-2 day-1
ds$GPP_mul_sqrt_VPD <- ds$GPP*sqrt(ds$VPD_F)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
ds$ET <- LE.to.ET(ds$LE_F_MDS,ds$TA_F)
ds$ET <- ds$ET * 24*3600  # mm s-1 to mm day-1 or kgH2O m-2 day-1
ds$Date <- date(ds$DateTime)
#-----------------------------------
ds$EBR <- (ds$LE_F_MDS + ds$H_F_MDS) / (ds$NETRAD - ds$G_F_MDS)  # Energy Balance Ratio
#-------------------------------------------------------------------------------
ds$ET_corr <- NA_real_  
ds$uWUEp_pred   <- NA_real_
ds$uWUEp_ref_yr <- NA_real_  
years <- sort(unique(ds$Year))
uwuep_site_yr <- list()
for (yr in years) {
  ds_year <- ds[ds$Year == yr, ]
  # ------------------------------
  # Reference uWUEp for 0.9–1.1 EBR
  # ------------------------------
  ref_mask <- with(ds_year, uWUEp_Mask & !is.na(EBR) & EBR >= 0.9 & EBR <= 1.1)
  temp_ref <- data.frame(
    x = ds_year$ET[ref_mask],
    y = ds_year$GPP_mul_sqrt_VPD[ref_mask]
  )
  temp_ref <- temp_ref[complete.cases(temp_ref), ]
  if (nrow(temp_ref) >= 10) {
    qt_fit_ref <- quantile_reg(temp_ref$x, temp_ref$y, PolyDeg = 0, tau = 0.95, weights = NA)
    uWUEp_ref <- qt_fit_ref[[1]]
  } else {
    uWUEp_ref <- NA_real_
  }
  ds$uWUEp_ref_yr[ds$Year == yr] <- uWUEp_ref
  
  # ------------------------------
  # Compute local uWUEp vs EBR bins
  # ------------------------------
  valid_mask <- with(ds_year, !is.na(EBR))
  temp_all <- data.frame(
    EBR = ds_year$EBR[valid_mask],
    x = ds_year$ET[valid_mask],
    y = ds_year$GPP_mul_sqrt_VPD[valid_mask]
  )
  temp_all <- temp_all[complete.cases(temp_all), ]
  temp_all$ebr_bins <- cut(temp_all$EBR, breaks = seq(0, 1.2, by = 0.05), include.lowest = TRUE)
  uwuep_by_bin <- temp_all %>%
    group_by(ebr_bins) %>%
    summarise(
      EBR_center = mean(EBR, na.rm = TRUE),
      uWUEp = {
        sub_ds <- data.frame(x, y)
        sub_ds <- sub_ds[complete.cases(sub_ds), ]
        if (nrow(sub_ds) >= 10) {
          qt_fit <- quantile_reg(sub_ds$x, sub_ds$y, PolyDeg = 0, tau = 0.95, weights = NA)
          qt_fit[[1]]
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    )
  
  uwuep_by_bin <- uwuep_by_bin[!is.na(uwuep_by_bin$uWUEp), ]
  if (nrow(uwuep_by_bin) > 0) {
    uwuep_by_bin$Year      <- yr
    uwuep_by_bin$uWUEp_ref <- uWUEp_ref
    uwuep_site_yr[[as.character(yr)]] <- uwuep_by_bin
  }
  # ------------------------------
  # LOESS smoothing
  # ------------------------------
  uwuep_by_bin <- uwuep_by_bin[complete.cases(uwuep_by_bin),]
  if (nrow(uwuep_by_bin) >= 5 & !is.na(uWUEp_ref)) {
    loess_fit <- loess(uWUEp ~ EBR_center, data = uwuep_by_bin, span = 0.75)
    # Predict uWUEp for each observation
    ds$uWUEp_pred[ds$Year == yr] <- predict(loess_fit, newdata = data.frame(EBR_center = ds$EBR[ds$Year == yr]))
    
    # Apply correction (uWUEp_pred / uWUEp_ref)
    ds$ET_corr[ds$Year == yr] <- with(ds[ds$Year == yr, ],
                                            ifelse(!is.na(uWUEp_pred) & uWUEp_pred > 0 & uWUEp_ref > 0 & (uWUEp_pred / uWUEp_ref >= 1 & uWUEp_pred / uWUEp_ref <= 5),
                                                   ET * (uWUEp_pred / uWUEp_ref),
                                                   NA_real_))
    ds$ET_corr[ds$Year == yr & ref_mask] <- ds$ET[ds$Year == yr & ref_mask]
  }
}
ds$LE_corr <- ET.to.LE(ds$ET_corr / (24*3600), ds$TA_F)
uwuep_site_yr <- bind_rows(uwuep_site_yr)
# ------------------------------------------------------------------------------
# Plot (LE original vs corrected LE)
plot(ds$LE_F_MDS, ds$LE_corr,
     xlab = expression("Original LE (W" ~ m^{-2} * ")"),
     ylab = expression("Corrected LE (W" ~ m^{-2} * ")"),
     pch = 16, col = "darkgray")
abline(0, 1, col = "red", lwd = 2)

