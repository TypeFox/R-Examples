# loading functions for input-oriented DEA

### RTS test statistic (4.6) in Simar and Wilson (2002)
RTSStatistic46 <- function(dist_crs, dist_vrs) {
  return( sum(dist_crs) / sum(dist_vrs) )
}

### RTS test statistic (4.5) in Simar and Wilson (2002)
RTSStatistic45 <- function(dist_crs, dist_vrs) {
  return( mean(dist_crs/dist_vrs) )
}

### RTS statistic (48) in Simar and Wilson (2011) J Prod Anal p.47
RTSStatistic48 <- function(dist_crs, dist_vrs) {
  return( mean(dist_crs/dist_vrs)-1 )
}

fixaround1 <- function(x) {
  x[ 1 - 1e-5 < x & x < 1 + 1e-5] <- 1
  x
}


######## Return to Scale Test (Simar and Wilson, 2002) ########

###### [ Inputs to the method ] ######
### X[firm, xfeat] - production inputs
### Y[firm, yfeat] - production output
### W[firm, xfeat] - input prices,  only necessary when model="costmin"
### model - either "input" (default) or "output"
### H0    - which RTS test to perform: "constant" vs variable (default), "non-increasing" vs variable.
### bw - theta bandwidth calculation method.
###      Either user supplied function or:
###      1) "silverman", "bw.nrd0" for Silverman rule
###      2) "rule" for rule of thumb
###      3) "cv", "bw.ucv" for unbiased cross-validation
###      Default is CV. 
### B     - number of bootstraps, default 2000
### alpha - confidence level, i.e., 0.05 (default)

###### [ Outputs from the method ] ######
### out$theta_H0_hat - theta, estimated reciprocal of DEA distance under CRS
### out$theta_vrs_hat - theta, estimated reciprocal of DEA distance under VRS
### out$w_hat      - estimated test statistic (4.6)
### out$w45_hat    - estimated test statistic (4.5), not used in the function, just calculated
### out$w_hat_boot - bootstrapped values for test statistic (4.6)
### out$w45_hat_boot - bootstrapped values for test statistic (4.5)
### out$pvalue     - p-value (from bootstrap)
### out$H0level    - test statistic value at alpha significance for data being CRS
### out$H0reject   - TRUE if H0 is rejected, i.e., pvalue < alpha
rts.test <- function(X, Y, W=NULL, model, H0="constant", bw="cv", B=2000, alpha=0.05) {
  out = list()
  if (! is.matrix(X)) X = as.matrix(X)
  if (! is.matrix(Y)) Y = as.matrix(Y)
  if (! is.numeric(X)) stop("X has to be numeric.")
  if (! is.numeric(Y)) stop("Y has to be numeric.")
  if ( any(is.na(X)) ) stop("X contains NA. Missing values are not supported.")
  if ( any(is.na(Y)) ) stop("Y contains NA. Missing values are not supported.")
  
  if (1/alpha > B) stop(sprintf("The number of bootstraps B (%d) has to be bigger than 1/alpha (%1.0f) to allow for confidence interval estimation.", B, 1/alpha))
  
  if (model == "costmin") {
    if ( is.null(W) )    stop("For costmin input prices (W) are necessary.")
    if (! is.matrix(W)) W = as.matrix(W)
    if (! is.numeric(W)) stop("W has to be numeric.")
    if ( any(is.na(W)) ) stop("W contains NA. Missing values are not supported.")
  }
  
  if (missing(model) || ! model %in% c("input", "output", "costmin") ) {
    stop("model has to be either 'input', 'output' or 'costmin'.")
  }
  
  if (! H0 %in% c("constant", "non-increasing")) {
    stop("H0 has to be either 'constant' or 'non-increasing'.")
  }
  
  if (! bw %in% c("cv", "silverman")) {
    stop("bw has to be either 'cv' (cross-validation) or 'silverman'.")
  }
  
  N    = nrow( X )
  Ydim = ncol( Y )
  
  if (model == "input")   dea = dea.input.rescaling
  if (model == "output")  dea = dea.output.rescaling
  if (model == "costmin") dea = dea.input.rescaling
  ## TODO:
  ## 1) add costmin model here
  ## 2) costmin uses theta_H0_hat and theta_vrs_hat from input oriented model
  ## 3) step (2) inside bootstrap uses input-oriented results for rescaling
  ## 4) step (4) inside bootstrap uses costmin dea to find results

  # (1) calculating input-oriented DEA with CRS and VRS:
  theta_H0_hat = fixaround1( as.vector( dea(XREF=X, YREF=Y, X=X, Y=Y, RTS=H0) ) )
  theta_vrs_hat = fixaround1( as.vector( dea(XREF=X, YREF=Y, X=X, Y=Y, RTS="variable") ) )
  w_hat     = RTSStatistic46( theta_H0_hat, theta_vrs_hat )
  
  delta_crs_hat = 1 / theta_H0_hat
  var_delta_crs_hat = var(delta_crs_hat)
  mean_delta_crs_hat = mean(delta_crs_hat)
  
  out$w_hat   = w_hat
  #out$w45_hat = RTSStatistic45( theta_H0_hat, theta_vrs_hat )
  out$w48_hat = RTSStatistic48( theta_H0_hat, theta_vrs_hat )
  out$theta_H0_hat = theta_H0_hat
  out$theta_vrs_hat = theta_vrs_hat

  # finding the bandwidth for kernel sampling:
  if (is.function(bw)) {
    ## user-supplied function
    bw_value = bw( as.vector(delta_crs_hat) )
  } else if (bw %in% c("silverman", "bw.nrd0") ) {
    ## silverman (using iqr)
    bw_value = bw.nrd0( as.vector(delta_crs_hat) )
  } else if (bw %in% c("cv", "bw.ucv")) {
    ## unbiased cross-validated bw
    suppressWarnings({
      bw_value = bw.ucv( as.vector(delta_crs_hat) )
    })
  } else {
    stop( sprintf("Illegal bandwidth type '%s'.", bw) )
  }
  
  # bootstrap loop:
  w_hat_boot = matrix(0, B, 1)
  w45_hat_boot = matrix(0, B, 1)
  w48_hat_boot = matrix(0, B, 1)
  for (i in 1:B) {
    # (2) reflection method, sampling reciprocals of distance from smooth kernel function:
    theta_crs_boot = 1 / sampling_delta_with_reflection( N, delta_crs_hat, bw_value, var_delta_crs_hat, mean_delta_crs_hat )
    # (3) new output based on the efficiencies:
    #X_boot     = X * (theta_H0_hat / theta_crs_boot)
    rescaling = (theta_H0_hat / theta_crs_boot)
    # (4) DEA efficiencies for scores:
    if (model == "costmin") {
      theta_H0_hat_boot  = dea.costmin.rescaling(XREF=X, YREF=Y, X=X, Y=Y, W=W, rescaling=rescaling, RTS=H0)
      theta_vrs_hat_boot = dea.costmin.rescaling(XREF=X, YREF=Y, X=X, Y=Y, W=W, rescaling=rescaling, RTS="variable")
    } else {
      theta_H0_hat_boot  = dea(XREF=X, YREF=Y, X=X, Y=Y, rescaling=rescaling, RTS=H0)
      theta_vrs_hat_boot = dea(XREF=X, YREF=Y, X=X, Y=Y, rescaling=rescaling, RTS="variable")
    }
    w_hat_boot[i]   = RTSStatistic46( theta_H0_hat_boot, theta_vrs_hat_boot )
    w45_hat_boot[i] = RTSStatistic45( theta_H0_hat_boot, theta_vrs_hat_boot )
    w48_hat_boot[i] = RTSStatistic48( theta_H0_hat_boot, theta_vrs_hat_boot )
  }
  out$w_hat_boot   = w_hat_boot
  #out$w45_hat_boot = w45_hat_boot
  out$w48_hat_boot = w48_hat_boot
  
  # number of bootstrapped statistic values smaller than w_hat:
  w_hat_ind = sum(w_hat_boot < w_hat)
  # p-value, +1 is added to avoid 0 p-value (minimum value is 1/B)
  out$pvalue   = (w_hat_ind + 1) / B
  #out$pvalue45 = (sum(w45_hat_boot < out$w45_hat) + 1) / B
  out$pvalue48 = (sum(w48_hat_boot < out$w48_hat) + 1) / B # need to check, Simar and Wilson (2011) p.47
  # reject H0 (i.e., data is CRS) when alpha is bigger than p-value
  out$H0reject   = out$pvalue < alpha
  #out$H0reject45 = out$pvalue45 < alpha
  out$H0reject48 = out$pvalue48 < alpha
  # sorting and finding out the cut off statistic value for conficence alpha:
  out$H0level   = sort(w_hat_boot)[   floor(alpha*B) ]
  #out$H0level45 = sort(w45_hat_boot)[ floor(alpha*B) ]
  out$H0level48 = sort(w48_hat_boot)[ floor(alpha*B) ]
  # store H0:
  out$H0       = H0
  # bw:
  out$bw       = bw
  out$bw_value = bw_value
  
  return(out)
}
