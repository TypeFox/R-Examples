#' Robust Regression using One-Sided Huber Function
#'
#' This function performs robust regression using M-estimation
#' using the one-sided Huber function, with residuals truncated at Q / (data$gregwt-1)
#' where data$gregwt is the generalized regression weight.
#'
#' @param formula The regression formula (e.g. income ~ employment + old.turnover if income is survey variable and employment and old.turnover are auxiliary variables).
#' @param data A data frame including the variables in formula, and gregwt (generalized regression estimator weight), and regwt (weight to be used in regression - will be set to 1 if missing).
#' @param Q The tuning parameter where large Q corresponds to no outlier treatment, and small Q corresponds to many outliers being flagged.
#' @param Qname Gives a variable name on data which contains a separate tuning parameter Q for every observation (either Q or Qname should be specified but not both).
#' @param maxit The maximum number of iterations.
#' @param stop Set to T to open a browser window (for debugging purposes)
#' @details
#' Uses iteratively reweighted least squares.
#' @return The final linear model fit (an object of class "lm").
#' @references
#' Clark, R. G. (1995), "Winsorisation methods in sample surveys," Masters thesis, Australian National University, http://hdl.handle.net/10440/1031.
#'
#' Kokic, P. and Bell, P. (1994), "Optimal winsorizing cutoffs for a stratified finite population estimator," J. Off. Stat., 10, 419-435.
#' @examples
#' robust.lm.onesided(formula=y~x1+x2,data=survdat.example,Q=250)

robust.lm.onesided <- function( formula , data , Q , Qname , maxit=100,stop=F){
  if(stop) browser()
  irls.w <- NA
  y <- data[,as.character(formula)[2]]
  if(!missing(Qname)) Q <- data[,Qname]
  if(!("regwt" %in% names(data))) data$regwt <- 1
  robust.resid <- (y-median(y))/(1.483*(quantile(y,probs=0.75)-quantile(y,probs=0.25)))
  #data$irls.w <- pmin(robust.resid,2) / robust.resid * data$regwt
  data$irls.w <- data$regwt
  data$irls.w[robust.resid>2] <- 2 / robust.resid[robust.resid>2] * data$regwt[robust.resid>2]
  old.beta <- rep(0,(length(as.character(formula))-2))
  delta <- 1
  iter <- 1
  while((iter<=100)&(delta>=1e-6)){
    wfit <- lm(formula=formula,data=data,weights=irls.w)
    data$irls.w <- pmin( y - wfit$fitted.values , Q / (data$gregwt-1) ) / ( y - wfit$fitted.values ) * data$regwt
    delta <- sqrt(sum((wfit$coef-old.beta)^2))
    old.beta <- wfit$coef
    iter <- iter+1
  }
  wfit
}
