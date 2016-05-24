####################
## Filter methods ##
####################

#' @aliases SR sMC LW RC
#' @title Filter methods for variable selection with Partial Least Squares.
#'
#' @description Various filter methods extracting and using information from 
#' \code{mvr} objects to assign importance to all included variables. Available 
#' methods are Significance Multivariate Correlation (sMC), Selectivity Ratio (SR), 
#' Variable Importance in Projections (VIP), Loading Weights (LW), Regression Coefficients (RC).
#'
#' @param pls.object \code{mvr} object from PLS regression.
#' @param opt.comp optimal number of components of PLS model.
#' @param p number of variables in PLS model.
#' @param X data matrix used as predictors in PLS modelling.
#' @param alpha_mc quantile significance for automatic selection of variables in \code{sMC}.
#'
#' @return A vector having the same lenght as the number of variables in the associated
#' PLS model. High values are associated with high importance, explained variance or
#' relevance to the model.
#'
#' @author Tahir Mehmood, Kristian Hovde Liland, Solve Sæbø.
#'
#' @references T. Mehmood, K.H. Liland, L. Snipen, S. Sæbø, A review of variable selection 
#' methods in Partial Least Squares Regression, Chemometrics and Intelligent Laboratory Systems
#' 118 (2012) 62-69.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{spa_pls}}, 
#' \code{\link{stpls}}, \code{\link{truncation}}, \code{\link{bve_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{ipw_pls}}, \code{\link{ga_pls}}, \code{\link{rep_pls}}.
#'
#' @examples
#' data(gasoline, package = "pls")
#' library(pls)
#' pls  <- plsr(octane ~ NIR, ncomp = 10, validation = "LOO", data = gasoline)
#' comp <- which.min(pls$validation$PRESS)
#' X    <- gasoline$NIR
#' vip <- VIP(pls, comp)
#' sr  <- SR (pls, comp, X)
#' smc <- sMC(pls, comp, X)
#' lw  <- LW (pls, comp)
#' rc  <- RC (pls, comp)
#' matplot(scale(cbind(vip, sr, smc, lw, rc)), type = 'l')
#'
#' @export
VIP <- function(pls.object, opt.comp, p = dim(pls.object$coef)[1]){
  # Variable importance in prediction
  W  <- pls.object$loading.weights
  Q  <- pls.object$Yloadings
  TT <- pls.object$scores
  Q2 <- as.numeric(Q)* as.numeric(Q)
  WW <- W*W/apply(W,2, function(x) sum(x*x))
  
  VIP <- sqrt(p*apply(as.matrix(Q2[1:opt.comp]*diag(crossprod(TT))[1:opt.comp]*WW[,1:opt.comp]),1,sum)/sum(Q2[1:opt.comp]*diag(crossprod(TT))[1:opt.comp]))
  VIP
}

#' @rdname VIP
#' @export
SR <- function(pls.object, opt.comp, X){
  # Selectivity ratio
  X   <- as.matrix(X)
  RC  <- pls.object$coefficients[,1,opt.comp]
  Wtp <- RC/norm(matrix(RC))
  Ttp <- X%*%Wtp
  Ptp <- (t(X)%*%Ttp)/c((t(Ttp)%*%Ttp))
  
  Xtp <- Ttp%*%t(Ptp)
  Xr  <- X-Xtp
  SR  <- colSums(Xtp*Xtp)/colSums(Xr*Xr)
  #  var.test(Xtp, Xr)
  SR
}

#' @rdname VIP
#' @export
sMC <- function(pls.object, opt.comp, X, alpha_mc = 0.05){
  # Significance Multivariate Correlation
  # [smcF smcFcrit SSCregression SSResidual] = smc(b, X)
  # Output:
  # smcF : SMC F-values for the list of variables
  # smcFcrit: F-critical cutoff threshold value for significant important variables (smcF>smcFcrit)
  #
  # In case of publication of any application of this method,
  # please, cite the original work:
  # T.N. Tran*, N.L. Afanador, L.M.C. Buydens, L. Blanchet, 
  # Interpretation of variable importance in Partial Least Squares with Significance Multivariate Correlation (sMC), 
  # Chemometrics and Intelligent Laboratory Systems, Volume 138, 15 November 2014, Pages 153-160
  # DOI: http://dx.doi.org/10.1016/j.chemolab.2014.08.005
  
  b  <- pls.object$coefficients[,1,opt.comp]
  X   <- as.matrix(X)
  
  n <- dim(X)[1]
  
  yhat <- X%*%b
  Xhat <- tcrossprod(yhat,b)/crossprod(b)[1]
  Xresidual <- X - Xhat
  
  SSCregression <- colSums(Xhat^2)
  SSResidual    <- colSums(Xresidual^2)
  
  MSCregression <- SSCregression # 1 degrees of freedom
  MSResidual    <- SSResidual/(n-2)
  
  smcF     <- MSCregression/MSResidual;
  smcFcrit <- qf(1-alpha_mc,1,n-2)
#  list(smcF=smcF, smcFcrit=smcFcrit)
  attr(smcF, "quantile") <- smcFcrit
  smcF
}

#' @rdname VIP
#' @export
LW <- function(pls.object, opt.comp)
  # Loading weights
  pls.object$loading.weights[ , opt.comp]

#' @rdname VIP
#' @export
RC <- function(pls.object, opt.comp)
  # Regression coefficients
  pls.object$coef[ , 1, opt.comp]

# Remove names and sort vector
simplify <- function(X){
  names(X) <- NULL
  sort(X)
}


###############
## filterPLS ##
###############

#' @title Optimisation of filters for Partial Least Squares
#'
#' @description Extract the index of influential variables based on threshold defiend for
#' LW (loading weights), RC (regression coef), JT (jackknife testing) and VIP (variable 
#' importance on projection).
#'
#' @param y vector of response values (\code{numeric} or \code{factor}).
#' @param X numeric predictor \code{matrix}.
#' @param ncomp integer number of components (default = 10).
#' @param ncomp.opt use the number of components corresponding to minimum error (minimum)
#'  or \code{ncomp} (same).
#' @param validation type of validation in the PLS modelling (default = "LOO").
#' @param LW.threshold threshold for Loading Weights if applied (default = NULL).
#' @param RC.threshold threshold for Regression Coefficients if applied (default = NULL).
#' @param JT.threshold threshold for Jackknife Testing if applied (default = NULL).
#' @param VIP.threshold threshold for Variable Importance on Projections if applied (default = NULL).
#' @param SR.threshold threshold for Selectivity Ration if applied (default = NULL).
#' @param sMC.threshold threshold for Significance Multivariate Correlation if applied (default = NULL).
#' @param ... additional paramters for \code{pls}, e.g. segmentation or similar.
#'
#' @details Filter methods are applied for variable selection with PLSR. This function can 
#' return selected variables and Root Mean Squared Error of Cross-Validation for various 
#' filter methods and determine optimum numbers of components.
#' 
#' @return Returns a list of lists containing filters (outer list), their selected variables,
#' optimal numbers of components and prediction accuracies.
#'
#' @author Tahir Mehmood, Kristian Hovde Liland, Solve Sæbø.
#'
#' @references T. Mehmood, K.H. Liland, L. Snipen, S. Sæbø, A review of variable selection 
#' methods in Partial Least Squares Regression, Chemometrics and Intelligent Laboratory Systems
#' 118 (2012) 62-69.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{spa_pls}}, 
#' \code{\link{stpls}}, \code{\link{truncation}}, \code{\link{bve_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{ipw_pls}}, \code{\link{ga_pls}}, \code{\link{rep_pls}}.
#'
#' @examples
#' data(gasoline, package = "pls")
#' \dontrun{
#' with( gasoline, filterPLSR(octane, NIR, ncomp = 10, "minimum", validation = "LOO",
#'  RC.threshold = c(0.1,0.5), SR.threshold = 0.5))
#' }
#' 
#' @import pls
#' @export
filterPLSR <- function(y, X, ncomp = ncomp, ncomp.opt = c("minimum","same"), validation = "LOO", LW.threshold = NULL, RC.threshold = NULL, JT.threshold = NULL, VIP.threshold = NULL, SR.threshold = NULL, sMC.threshold = NULL,...){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # Primary PLSR fitting
  if(!is.null(JT.threshold)){
    pls.object <- plsr(y ~ X, ncomp = ncomp, validation = validation, jackknife = TRUE, ...)
  } else {
    pls.object <- plsr(y ~ X, ncomp = ncomp, validation = validation, ...)
  }
  Press <- pls.object$valid$PRESS[1,]
  opt.comp <- which.min(Press)
  
  # Apply filter methods
  selections <- list()
  if(!is.null(LW.threshold)){
    LWvalues <- LW(pls.object,opt.comp)
    if(is.logical(LW.threshold))
      LW.threshold <- 0.05
    if(length(LW.threshold) == 1){ # Single threshold
      selections$LW <- simplify(which(LWvalues > LW.threshold))
      if(length(selections$LW) < ncomp)
        selections$LW <- simplify(order(LWvalues, decreasing = TRUE)[1:ncomp])
    } else { # Optimise over several thresholds
      selections$LW <- list()
      selections$LW$comps  <- numeric(length(LW.threshold))
      selections$LW$RMSECV <- numeric(length(LW.threshold))
      names(selections$LW$comps)  <- LW.threshold
      names(selections$LW$RMSECV) <- LW.threshold
      for(i in 1:length(LW.threshold)){
        selections$LW[[i+2]] <- simplify(which(LWvalues > LW.threshold[i]))
        if(length(selections$LW[[i+2]]) < ncomp)
          selections$LW[[i+2]] <- simplify(order(LWvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$LW[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
        selections$LW$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
        selections$LW$comps[i]  <- comp
      }
    }
  }
  
  if(!is.null(RC.threshold)){
    RCvalues <-  RC(pls.object,opt.comp)
    if(is.logical(RC.threshold))
      RC.threshold <- 0.01
    if(length(RC.threshold) == 1){ # Single threshold
      selections$RC <- simplify(which(RCvalues > RC.threshold))
      if(length(selections$RC) < ncomp)
        selections$RC <- simplify(order(RCvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$RC <- list()
      selections$RC$comps  <- numeric(length(RC.threshold))
      selections$RC$RMSECV <- numeric(length(RC.threshold))
      names(selections$RC$comps)  <- RC.threshold
      names(selections$RC$RMSECV) <- RC.threshold
      for(i in 1:length(RC.threshold)){
        selections$RC[[i+2]] <- simplify(which(RCvalues > RC.threshold[i]))
        if(length(selections$RC[[i+2]]) < ncomp)
          selections$RC[[i+2]] <- simplify(order(LWvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$RC[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
        selections$RC$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
        selections$RC$comps[i]  <- comp
      }
    }
  }
  
  if(!is.null(VIP.threshold)){
    VIPvalues <- VIP(pls.object, opt.comp)
    if(is.logical(VIP.threshold))
      VIP.threshold <- 1
    if(length(VIP.threshold) == 1){ # Single threshold
      selections$VIP <- simplify(which(VIPvalues > VIP.threshold))
      if(length(selections$VIP) < ncomp)
        selections$VIP <- simplify(order(VIPvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$VIP <- list()
      selections$VIP$comps  <- numeric(length(VIP.threshold))
      selections$VIP$RMSECV <- numeric(length(VIP.threshold))
      names(selections$VIP$comps)  <- VIP.threshold
      names(selections$VIP$RMSECV) <- VIP.threshold
      for(i in 1:length(VIP.threshold)){
        selections$VIP[[i+2]] <- simplify(which(VIPvalues > VIP.threshold[i]))
        if(length(selections$VIP[[i+2]]) < ncomp)
          selections$VIP[[i+2]] <- simplify(order(VIPvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$VIP[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
        selections$VIP$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
        selections$VIP$comps[i]  <- comp
      }
    }
  }
  
  if(!is.null(SR.threshold)){
    SRvalues <- SR(pls.object, opt.comp, X)
    if(is.logical(SR.threshold))
      SR.threshold <- pf(0.99, df1 = n-1, df2 = n-2)
    if(length(SR.threshold) == 1){ # Single threshold
      selections$SR <- simplify(which(SRvalues > SR.threshold))
      if(length(selections$SR) < ncomp)
        selections$SR <- simplify(order(SRvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$SR <- list()
      selections$SR$comps  <- numeric(length(SR.threshold))
      selections$SR$RMSECV <- numeric(length(SR.threshold))
      names(selections$SR$comps)  <- SR.threshold
      names(selections$SR$RMSECV) <- SR.threshold
      for(i in 1:length(SR.threshold)){
        selections$SR[[i+2]] <- simplify(which(SRvalues > SR.threshold[i]))
        if(length(selections$SR[[i+2]]) < ncomp)
          selections$SR[[i+2]] <- simplify(order(LWvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$SR[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
        selections$SR$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
        selectinos$SR$comps[i]  <- comp
      }
    }
  }
  
  if(!is.null(sMC.threshold)){
    sMCvalues <- sMC(pls.object, opt.comp, X)$smcF
    if(is.logical(sMC.threshold))
      sMC.threshold <- pf(0.99, df1 = n-1, df2 = n-2)
    if(length(sMC.threshold) == 1){ # Single threshold
      selections$sMC <- simplify(which(sMCvalues > sMC.threshold))
      if(length(selections$sMC) < ncomp)
        selections$sMC <- simplify(order(sMCvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$sMC <- list()
      selections$sMC$comps  <- numeric(length(sMC.threshold))
      selections$sMC$RMSECV <- numeric(length(sMC.threshold))
      names(selections$sMC$comps)  <- sMC.threshold
      names(selections$sMC$RMSECV) <- sMC.threshold
      for(i in 1:length(sMC.threshold)){
        selections$sMC[[i+2]] <- simplify(which(sMCvalues > sMC.threshold[i]))
        if(length(selections$sMC[[i+2]]) < ncomp)
          selections$sMC[[i+2]] <- simplify(order(LWvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$sMC[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
        selections$sMC$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
        selections$sMC$comps[i]  <- comp
      }
    }
  }
  
  if(!is.null(JT.threshold)){
    JTvalues <- jack.test(pls.object, ncomp = opt.comp)$pvalues[,1,1]
    if(is.logical(JT.threshold))
      JT.threshold <- 0.05
    if(length(JT.threshold) == 1){ # Single threshold
      selections$JT <- simplify(which(JTvalues > JT.threshold))
      if(length(selections$JT) < ncomp)
        selections$JT <- simplify(order(JTvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$JT <- list()
      selections$JT$comps  <- numeric(length(JT.threshold))
      selections$JT$RMSECV <- numeric(length(JT.threshold))
      names(selections$JT$comps)  <- JT.threshold
      names(selections$JT$RMSECV) <- JT.threshold
      for(i in 1:length(JT.threshold)){
        selections$JT[[i+2]] <- simplify(which(JTvalues > JT.threshold[i]))
        if(length(selections$JT[[i+2]]) < ncomp)
          selections$JT[[i+2]] <- simplify(order(LWvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$JT[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
        selections$JT$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
        selections$JT$comps[i]  <- comp
      }
    }
  }
  
  selections
}
