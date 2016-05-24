# this file does bootstrap by Simar and Wilson (2007)
############### Algorithm 2 (Page 42) ###############

### Inputs to the method ###
### Xbar[firm, xfeat] - production inputs multiplied by prices
### Y[firm, yfeat] - production outputs
### Z[firm, zfeat] - environmental variables
### RTS       - returns to scale: "variable" (default), "non-increasing", "constant"
### L1/L2     - number of bootstraps
### alpha     - confidence interval, i.e., 0.05
#costmin.tone.bias.correction.sw07 <- function(Xbar, Y, Z, RTS="variable", L1=100, L2=1000, alpha=0.05) {
#  ## calling the general bias correction method:
#  bias.correction.sw07(X=Xbar, Y=Y, Z=Z, RTS=RTS, L1=L1, L2=L2, alpha=alpha, 
#                       deaMethod=costmin.aggr.rglpk.scaling)
#}


### Inputs to the method ###
### X[firm, xfeat] - production inputs in units
### W[firm, xfeat] - prices of inputs for each firm
### Y[firm, yfeat] - production outputs
### Z[firm, zfeat] - environmental variables
### RTS       - returns to scale: "variable" (default), "non-increasing", "constant"
### L1/L2     - number of bootstraps
### alpha     - confidence interval, i.e., 0.05
dea.env.robust <- function(X, Y, W=NULL, Z, model, RTS="variable", 
                           L1=100, L2=2000, alpha=0.05) {
  if (missing(model) || ! model %in% c("input", "output", "costmin", "costmin-tone") ) {
    stop("Parameter 'model' has to be either 'input', 'output', 'costmin' or 'costmin-tone'.")
  }
  if (1/alpha > L2) stop(sprintf("The number of bootstraps L2 (%d) has to be bigger than 1/alpha (%1.0f) to allow for confidence interval estimation.", L2, 1/alpha))
  
  if ( ! is.matrix(X)) { X = as.matrix(X) }
  if ( ! is.matrix(Y)) { Y = as.matrix(Y) }
  if ( ! is.matrix(Z)) { Z = as.matrix(Z) }
  if ( ! is.numeric(X)) { stop("X has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(Y)) { stop("Y has to be numeric matrix or data.frame.") }
  if ( ! is.numeric(Z)) { stop("Z has to be numeric matrix or data.frame.") }
  if ( any(is.na(X)) ) stop("X contains NA. Missing values are not supported.")
  if ( any(is.na(Y)) ) stop("Y contains NA. Missing values are not supported.")
  if ( any(is.na(Z)) ) stop("Z contains NA. Missing values are not supported.")
  if (nrow(X) != nrow(Y)) { stop( sprintf("Number of rows in X (%d) does not equal the number of rows in Y (%d)", nrow(X), nrow(Y)) ) }
  if (nrow(X) != nrow(Z)) { stop( sprintf("Number of rows in X (%d) does not equal the number of rows in Z (%d)", nrow(X), nrow(Z)) ) }
  
  if (model == "input") {
    return( bias.correction.sw07(X=X, Y=Y, Z=Z, RTS=RTS, L1=L1, L2=L2, alpha=alpha, 
                                 deaMethod=dea.input.rescaling))
  } else if (model == "output") {
    return( bias.correction.sw07(X=X, Y=Y, Z=Z, RTS=RTS, L1=L1, L2=L2, alpha=alpha, 
                                 deaMethod=dea.output.rescaling))
  } else {
    ## costmin (Fare's costmin)
    if ( is.null(W)) { stop("W (input prices) has to be numeric matrix or data.frame.") }
    if ( ! is.matrix(W)) { W = as.matrix(W) }
    if ( ! is.numeric(W)) { stop("W (input prices) has to be numeric matrix or data.frame.") }
    if ( any(is.na(W)) ) stop("W (input prices) contains NA. Missing values are not supported.")
    if (nrow(X) != nrow(W)) { stop( sprintf("Number of rows in X (%d) does not equal the number of rows in W (%d)", nrow(X), nrow(W)) ) }
    if (ncol(X) != ncol(W)) { stop( sprintf("Number of columns in X (%d) does not equal the number of columns in W (%d)", ncol(X), ncol(W)) ) }

    ## TODO: add if for costmine-tone here
    if (model == "costmin-tone") stop( "costmin-tone is not available yet." )
    
    dm <- function(XREF, YREF, X, Y, RTS, rescaling=1) {
      dea.costmin.rescaling(XREF=XREF, YREF=YREF, W=W, X=X, Y=Y, RTS=RTS, rescaling=rescaling)
    }
    ## calling the general bias correction method:
    return( bias.correction.sw07(X=X, Y=Y, Z=Z, RTS=RTS, L1=L1, L2=L2, alpha=alpha, deaMethod=dm) )
    
  }
}


### Inputs to the method ###
### X[firm, xfeat] - production inputs
### Y[firm, yfeat] - production outputs
### Z[firm, zfeat] - environmental variables
### RTS       - returns to scale: "variable", "non-increasing", "constant"
### L1/L2     - number of bootstraps
### alpha     - confidence interval, i.e., 0.05
### deaMethod - method to use for calculating efficiency
bias.correction.sw07 <- function(X, Y, Z, RTS, L1, L2, alpha, deaMethod) {
  out = list()
  
  if (!is.matrix(X)) { X = as.matrix(X) }
  if (!is.matrix(Y)) { Y = as.matrix(Y) }
  if (!is.matrix(Z)) { Z = as.matrix(Z) }
  
  # add constant for easy matrix multiplication in truncated regression
  Z1 = cbind(1.0, Z)
  
  # number of firms
  N    = nrow( as.matrix(X) )
  
  # 1. Compute the original DEA (input oriented)
  # now we get delta_hat
  delta_hat = 1.0 / deaMethod(XREF=X, YREF=Y, X=X, Y=Y, RTS=RTS)
  # removing numerical errors
  delta_hat[ 1 - 1e-5 < delta_hat & delta_hat < 1 ] = 1
  out$delta_hat = delta_hat
  
  # 2. Maximum likelihood for estimating beta and sigma_{\epsilon}
  # - only use delta_hat>1
  # from regression we get beta_hat and sigma_hat
  
  ind_not1 = out$delta_hat>1
  ty   = out$delta_hat[ind_not1]
  tx   = Z[ind_not1,,drop=F]
  if (sum(ind_not1) <= ncol(tx) + 1) {
    ## too few non-frontier samples
    stop(paste0(
      "Too few firms to build the model for environmental variables.\n",
      sprintf("The number of firms not on the frontier is %d ", sum(ind_not1)),
      sprintf("and the number of environmental variables is %d + 1 constant term.", ncol(tx)),
      "\n\n",
      "To solve this issue either:",
      "\n1) Increase the number of firms (which should increase the number of non-frontier firms).",
      "\n2) Decrease the number of environmental variables.",
      "\n3) Decrease the number of output variables, which should move some firms off the frontier."
    ))
  }
  dataf = data.frame(ty = ty, tx)
  suppressWarnings({ mlYZ = truncreg::truncreg(ty ~ ., data=dataf, point=1.0, direction="left") })
  beta_hat      = mlYZ$coefficients[1:(length(mlYZ$coefficients)-1)]
  sigma_hat     = mlYZ$coefficients[length(mlYZ$coefficients)]
  out$beta_hat  = beta_hat
  out$sigma_hat = sigma_hat
  
  Zbeta = Z1 %*% matrix(beta_hat)
  
  # 3. Loop of L1 times (100) to get bootstrapped values of delta_hat for each firm
  delta_hat_star = matrix(0, N, L1)
  for (i in 1:L1) {
    # 3.1 Draw \epsilon_i for each observation from trunc. normal
    epsilon = rtruncnorm(N, a=1-Zbeta, sd=sigma_hat )
    
    # 3.2 Compute \delta^*_i
    # delta_star = t(Z)*beta_hat + epsilon
    delta_star = Zbeta + epsilon
    
    # 3.3 Compute y^*, rescaling
    # for i=1:n
    rescaleFactor = as.vector( delta_star/delta_hat )
    #X_star = X * ( delta_star/delta_hat )[,]
    
    # 3.4 Compute DEA for y^* and x to get \widehat{\delta}^*_i
    #dea.boot = input.dea.rglpk(XREF=X_star, YREF=Y, X=X, Y=Y, RTS=RTS)
    #delta_hat_star[,i] = 1.0 / dea.boot$thetaOpt
    delta_hat_star[,i] = 1.0 / deaMethod(XREF=X, YREF=Y, X=X, Y=Y, RTS=RTS, 
                                         rescaling=rescaleFactor)
  }
  
  #  delta_ci = apply( delta_hat_star, 1, quantile, c(alpha/2, 1-alpha/2) )
  #  out$delta_ci_low  = delta_ci[1,] - 2*bias
  #  out$delta_ci_high = delta_ci[2,] - 2*bias
  
  # 4. for each firm compute bias-corrected estimator \doublehat{\delta}_i.
  # delta_hat_hat = mean of bootstrap values for each firm
  # bias[i] = mean of delta_hat_star for firm i - t(Z[:,i])*beta_hat
  # delta_hat_hat = delta_hat - bias
  
  #bias = rowMeans(delta_hat_star) - Zbeta
  delta_hat_hat = 2*delta_hat - rowMeans(delta_hat_star)
  
  # list for output variables:
  out$bias          = rowMeans(delta_hat_star) - delta_hat
  out$beta_hat      = beta_hat
  out$delta_hat_hat = delta_hat_hat
  
  # estimate CI for delta
  # c(alpha/2, 1-alpha/2) is c(0.025, 0.975) in case alpha = 0.05.
  delta_ci = apply( delta_hat_star, 1, quantile, c(alpha/2, 1-alpha/2) )
  out$delta_ci_low  = delta_ci[1,] - 2*out$bias
  out$delta_ci_high = delta_ci[2,] - 2*out$bias

  ## from Kneip et al (2008), page 1676-1677:
  ## (there is theta, we have delta = 1 / theta)
  ratio_ci_low  = out$delta_hat / delta_ci[1,] - 1
  ratio_ci_high = out$delta_hat / delta_ci[2,] - 1

  out$delta_ci_kneip_low  = 1 / ((1/out$delta_hat) / (1 + ratio_ci_high))
  out$delta_ci_kneip_high = 1 / ((1/out$delta_hat) / (1 + ratio_ci_low))

  # 5. Maximum likelihood to estimate \doublehat{\beta} and \doublehat{\sigma}
  # - only use delta_hat>1
  # - Z should scaled by mean(Z)
  # from regression we get beta_hat_hat and sigma_hat_hat
  ind_not1 = delta_hat_hat>1
  out$delta_hat_hat_used = mean(ind_not1)
  
  ty   = delta_hat_hat[ind_not1]
  tx   = Z[ind_not1,]
  dataf= data.frame(ty = ty, tx)
  suppressWarnings({ mlYZ = truncreg::truncreg(ty ~ ., data=dataf, point=1.0, direction="left") })
  beta_hat_hat     = mlYZ$coefficients[ 1:(length(mlYZ$coefficients)-1) ]
  sigma_hat_hat    = mlYZ$coefficients[ length(mlYZ$coefficients) ]
  out$beta_hat_hat = beta_hat_hat
  out$sigma_hat_hat= sigma_hat_hat
  
  # 6. Bootstrap loop for L2 (2000) steps.
  beta_hat_hat_boot  = matrix(0, length(beta_hat), L2)
  sigma_hat_hat_boot = matrix(0, 1, L2)
  
  Zbeta_hat_hat = Z1 %*% matrix(beta_hat_hat)
  out$Zbeta_hat_hat = Zbeta_hat_hat
  
  for (i in 1:L2) {
    # 6.1 Draw \epsilon_i for each observation from trunc. normal
    epsilon = rtruncnorm(N, a=1-Zbeta_hat_hat, sd=sigma_hat_hat )
    
    # 6.2 Compute \delta^{**}_i
    # delta_star_star = t(Z)*beta_hat_hat + epsilon
    delta_star_star = Zbeta_hat_hat + epsilon
    
    # 6.3. Maximum likelihood to estimate \doublehat{\beta} and \doublehat{\sigma}
    # use FEAR::treg(Y=delta_star_star, X=Z)
    # - only use delta_star_star>1
    # - Z should scaled by mean(Z)
    # from regression we get beta_hat_hat and sigma_hat_hat
    ind_not1 = delta_star_star>1.0
    
    #mlYZ = FEAR::treg(Y=delta_star_star[ind_not1], X=Z[ind_not1,])
    #if (mlYZ$ier!=0) {
    #  stop("FEAR::treg could not solve the truncated regression (Step 6.3)!")
    #}
    #beta_hat_hat_star  = beta_hat_hat_star  + mlYZ$bhat
    #sigma_hat_hat_star = sigma_hat_hat_star + mlYZ$sighat
    
    # for package truncreg:
    ty = delta_star_star[ind_not1]
    tx = Z[ind_not1,]
    
    if (sum(is.infinite(ty))>0) {
      stop("Infinite values in step 6.3")
      #out$ty = data.frame(ty)
      #return(out)
    }
    dataf= data.frame(ty = ty, tx)
    suppressWarnings({ mlYZ = truncreg::truncreg(ty ~ ., data=dataf, point=1.0, direction="left") })
    beta_hat_hat_boot[,i]   = mlYZ$coefficients[1:(length(mlYZ$coefficients)-1)]
    sigma_hat_hat_boot[1,i] = mlYZ$coefficients[length(mlYZ$coefficients)]
    
  }
  # calculating the mean:
  beta_hat_hat_star = rowMeans( beta_hat_hat_boot )
  sigma_hat_hat_star = rowMeans( sigma_hat_hat_boot )
  
  out$beta_hat_hat_star  = beta_hat_hat_star
  out$sigma_hat_hat_star = sigma_hat_hat_star
  
  # 7. Calculate confidence interval for \doublehat{\beta}^* and \doublehat{\sigma}^*
  out$beta_ci  = t(apply( beta_hat_hat_boot, 1, quantile, c(alpha/2, 1-alpha/2) ))
  out$sigma_ci = t(apply( sigma_hat_hat_boot, 1, quantile, c(alpha/2, 1-alpha/2) ))
  
  ## centering CI to the original estimates \doublehat{\beta}^* and \doublehat{\sigma}^*
  out$beta_ci  = out$beta_ci  + as.matrix(beta_hat_hat - beta_hat_hat_star)[,]
  out$sigma_ci = out$sigma_ci + (sigma_hat_hat - sigma_hat_hat_star)
  
  rownames(out$beta_ci)  = rownames(out$beta_hat_hat)
  
  return(out)
}
