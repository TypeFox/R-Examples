### Bias-corrected DEA, input oriented case, (Simar and Wilson, 1998) 

###### [ Inputs to the method ] ######
### X[firm, input_feat]  - inputs
### Y[firm, output_feat] - outputs
### W[firm, input_feat]  - input prices, only needed for cost minimization
### RTS   - returns to scale, "variable" (default), "constant" and "non-increasing"
### B     - number of bootstraps (replications), default is 1000
### alpha - confidence level, default is 0.05

###### [ Outputs from the method ] ######
### out$theta_hat      - original DEA efficiency score
### out$theta_hat_star - efficiency scores for bootstrapped data
### out$bias           - bias
### out$theta_hat_hat  - bias corrected efficiency scores 
dea.robust <- function(X, Y, W=NULL, model, RTS="variable", B=1000, alpha=0.05,
                       bw = "bw.ucv", bw_mult = 1) {
  if (missing(model) || ! model %in% c("input", "output", "costmin") ) {
    stop("Parameter 'model' has to be either 'input', 'output' or 'costmin'.")
  }
  if (1/alpha > B) stop(sprintf("The number of bootstraps B (%d) has to be bigger than 1/alpha (%1.0f) to allow for confidence interval estimation.", B, 1/alpha))
  if ( ! is.matrix(X)) { X = as.matrix(X) }
  if ( ! is.matrix(Y)) { Y = as.matrix(Y) }
  if ( ! is.numeric(X)) { stop("X has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(Y)) { stop("Y has to be numeric matrix or data.frame.") }
  if ( any(is.na(X)) ) stop("X contains NA. Missing values are not supported.")
  if ( any(is.na(Y)) ) stop("Y contains NA. Missing values are not supported.")
  if ( bw_mult <= 0 )  stop("bw_mult has to be a positive value.")
    
  if (nrow(X) != nrow(Y)) { stop( sprintf("Number of rows in X (%d) does not equal the number of rows in Y (%d)", nrow(X), nrow(Y)) ) }
  
  if (model == "input") {
    return( bias.correction.sw98(X=X, Y=Y, RTS=RTS, B=B, alpha=alpha, bw=bw, bw_mult=bw_mult,
                                 deaMethod=dea.input.rescaling) )
  } else if (model == "output") {
    return( bias.correction.sw98(X=X, Y=Y, RTS=RTS, B=B, alpha=alpha, bw=bw, bw_mult=bw_mult,
                                 deaMethod=dea.output.rescaling) )
  } else { 
    return( dea.robust.costmin(X=X, Y=Y, W=W, RTS=RTS, B=B, alpha=alpha, bw=bw, bw_mult=bw_mult) )
  }
}

### inputs:
### bw - theta bandwidth calculation method, either user supplied function or
###      1) "silverman", "bw.nrd0" for Silverman rule
###      2) "rule" for rule of thumb
###      3) "cv", "bw.ucv" for unbiased cross-validation
bias.correction.sw98 <- function(X, Y, RTS, B, alpha, deaMethod, bw, bw_mult) {
  if (!is.matrix(X)) X = as.matrix(X)
  if (!is.matrix(Y)) Y = as.matrix(Y)
  out = list()

  # number of firms
  N    = nrow( as.matrix(X) )
  Ydim = ncol( as.matrix(Y) )
  Xdim = ncol( as.matrix(X) )
    
  # (1) calculating original theta_hat in DEA:
  theta_hat      = deaMethod(XREF=X, YREF=Y, X=X, Y=Y, RTS=RTS)
  delta_hat      = 1 / theta_hat
  var_delta_hat  = var(delta_hat)
  mean_delta_hat = mean(delta_hat)
  out$theta_hat  = theta_hat
  
  # finding the bandwidth for kernel sampling:
  delta_hat_m  = as.vector(delta_hat[delta_hat > 1])
  delta_hat_2m = c(delta_hat_m, 2 - delta_hat_m)
  if (is.function(bw)) {
    bw_value = bw( delta_hat_2m )
  } else if (bw == "rule") {
    bw_value = bandwidth_rule( ncol(X), ncol(Y), nrow(X) )
  } else if (bw == "rulesq") {
    bw_value = bandwidth_rule( ncol(X), ncol(Y), nrow(X) ) ^ 2
  } else if (bw %in% c("silverman", "bw.nrd0") ) {
    bw_value = bw.nrd0( delta_hat_2m )
  } else if (bw %in% c("cv", "bw.ucv")) {
    suppressWarnings({
      bw_value = bw.ucv( delta_hat_2m )
    })
  } else if (bw %in% c("cv.adjusted")) {
    ## Formula (4.86) from "The measurement of Productive Efficiency" (2008)
    suppressWarnings({
      bw_value = bw.ucv(delta_hat_2m)
      m        = length(delta_hat_m)
      n        = length(delta_hat)
      bw_value = bw_value * 2^(1/5) * (m / n)^(1/5) * sd(delta_hat) / sd(delta_hat_2m)
    })
  } else {
    stop( sprintf("Illegal bandwidth type '%s'.", bw) )
  }
  bw_value = bw_value * bw_mult
  
  # bootstrap loop (order of the loops is like in SW98 paper):
  theta_hat_star = matrix(0, N, B)
  for (b in 1:B) {
    # (2) reflection method, sampling reciprocals of distance from smooth kernel function:
    # use sampling from log normal to avoid negative values
    theta_hat_boot = 1 / sampling_delta_with_reflection( N, delta_hat, bw_value, var_delta_hat, mean_delta_hat )
    # (3) new output based on the efficiencies:
    rescaleFactor = as.vector(theta_hat / theta_hat_boot)
    
    # (4) DEA efficiencies for scores:
    # The efficiency of the original data under the bootstrapped technology:
    theta_hat_star[,b] = deaMethod(XREF=X, YREF=Y, X=X, Y=Y, RTS=RTS, 
                                   rescaling=rescaleFactor)
  }
  out$theta_hat_star = theta_hat_star
  out$bw = bw_value
  
  # (5) bias
  bias = rowMeans(theta_hat_star) - theta_hat
  # (6) bias-corrected estimator
  theta_hat_hat = theta_hat - bias

  out$bias          = bias
  out$theta_hat_hat = theta_hat_hat

  # calculating confidence interval for theta:
  theta_ci = apply( theta_hat_star, 1, quantile, c(alpha/2, 1-alpha/2) )
  out$theta_ci_low  = theta_ci[1,] - 2*bias
  out$theta_ci_high = theta_ci[2,] - 2*bias

  return(out)
}


### Inputs to the method ###
### Robust costmin for Fare's definition (Besstremyannaya et al, in press)
### X[firm, xfeat] - production inputs in units
### W[firm, xfeat] - prices of inputs for each firm
### Y[firm, yfeat] - production outputs
### RTS       - returns to scale: "variable" (default), "non-increasing", "constant"
### B         - number of bootstraps
### alpha     - confidence interval, i.e., 0.05
### bw        - "cv", "silverman", "rule"
### bw_mult   - bw multiplier, 1 means no effect
dea.robust.costmin <- function(X, Y, W, RTS, B, alpha, bw, bw_mult) {
  if (!is.matrix(X)) X = as.matrix(X)
  if (!is.matrix(Y)) Y = as.matrix(Y)
  if (!is.matrix(W)) W = as.matrix(W)
  if ( any(is.na(X)) ) stop("X contains NA. Missing values are not supported.")
  if ( any(is.na(Y)) ) stop("Y contains NA. Missing values are not supported.")
  if ( any(is.na(W)) ) stop("W contains NA. Missing values are not supported.")
  
  out = list()
  
  # number of firms
  N    = nrow( as.matrix(X) )
  Ydim = ncol( as.matrix(Y) )
  Xdim = ncol( as.matrix(X) )
  
  # (1) calculating original theta_hat in DEA:
  theta_hat     = dea.input.rescaling(XREF=X, YREF=Y, X=X, Y=Y, RTS=RTS, rescaling=1)
  gamma_hat     = dea.costmin.rescaling(XREF=X, YREF=Y, W=W, X=X, Y=Y, RTS=RTS, rescaling=1)
  out$theta_hat = theta_hat
  out$gamma_hat = gamma_hat
  delta_hat     = 1 / theta_hat
  var_delta_hat  = var(delta_hat)
  mean_delta_hat = mean(delta_hat)
  # finding the bandwidth for kernel sampling:
  delta_hat_m  = as.vector(delta_hat[delta_hat > 1])
  delta_hat_2m = c(delta_hat_m, 2 - delta_hat_m)
  if (is.function(bw)) {
    bw_value = bw( delta_hat_2m )
  } else if (bw == "rule") {
    bw_value = bandwidth_rule( ncol(X), ncol(Y), nrow(X) )
  } else if (bw == "rulesq") {
    bw_value = bandwidth_rule( ncol(X), ncol(Y), nrow(X) ) ^ 2
  } else if (bw %in% c("silverman", "bw.nrd0") ) {
    bw_value = bw.nrd0( delta_hat_2m )
  } else if (bw %in% c("cv", "bw.ucv")) {
    suppressWarnings({
      bw_value = bw.ucv( delta_hat_2m )
    })
  } else if (bw %in% c("cv.adjusted")) {
    ## Formula (4.86) from "The measurement of Productive Efficiency" (2008)
    suppressWarnings({
      bw_value = bw.ucv(delta_hat_2m)
      m        = length(delta_hat_m)
      n        = length(delta_hat)
      bw_value = bw_value * 2^(1/5) * (m / n)^(1/5) * sd(delta_hat) / sd(delta_hat_2m)
    })
  } else {
    stop( sprintf("Illegal bandwidth type '%s'.", bw) )
  }
  bw_value = bw_value * bw_mult
  
  out$bw = bw_value
  # bootstrap loop (order of the loops is like in SW98 paper):
  gamma_hat_star = matrix(0, N, B)
  #w_hat_boot     = matrix(0, B, 1)
  for (b in 1:B) {
    # (2) reflection method, sampling reciprocals of distance from smooth kernel function:
    theta_hat_boot = 1 / sampling_delta_with_reflection( N, delta_hat, bw_value, var_delta_hat, mean_delta_hat )
    # (3) new output based on the efficiencies:
    rescaleFactor = as.vector(theta_hat / theta_hat_boot)
    
    # (4) DEA efficiencies for scores:
    # The efficiency of the original data under the bootstrapped technology:
    gamma_hat_star[,b] = dea.costmin.rescaling(XREF=X, YREF=Y, W=W, X=X, Y=Y, RTS=RTS, 
                                               rescaling=rescaleFactor)
  }
  out$gamma_hat_star = gamma_hat_star
  
  # (5) bias
  bias = rowMeans(gamma_hat_star) - gamma_hat
  # (6) bias-corrected estimator 
  gamma_hat_hat = gamma_hat - bias
  
  out$bias          = bias
  out$gamma_hat_hat = gamma_hat_hat
  
  # calculating confidence interval for theta:
  gamma_ci = apply( gamma_hat_star, 1, quantile, c(alpha/2, 1-alpha/2) )
  out$gamma_ci_low  = gamma_ci[1,] - 2*bias
  out$gamma_ci_high = gamma_ci[2,] - 2*bias
  
  ## from Kneip et al (2008), page 1676-1677:
  ratio_ci_low  = gamma_ci[1,] / out$gamma_hat - 1
  ratio_ci_high = gamma_ci[2,] / out$gamma_hat - 1
  
  ## ratio ci low is used for gamma ci high
  out$gamma_ci_kneip_high = out$gamma_hat / (1 + ratio_ci_low) #sw98.fare
  out$gamma_ci_kneip_low  = out$gamma_hat / (1 + ratio_ci_high)
  
  ## paper approach
  #out$ratio2_ci = apply( gamma_hat_star / out$gamma_hat - 1, 1, quantile, c(alpha/2, 1-alpha/2) )
  
  
  return(out)
}
