arimapar <-
function(timeseries, na.action=na.omit, xreg=NULL){
  if(is.null(timeseries)) stop("timeseries is required and must have positive length")
  
  ts <- ts(na.action(timeseries))
  
  nobs <- length(ts)
  reg <- cbind(1:nobs,xreg)
  
  fit <- auto.arima(ts,xreg=ts(reg,start=1))
  
  #fit$arma -> A compact form of the specification, as a vector giving the number of AR, MA, seasonal AR and seasonal MA coefficients, plus the period and the number of non-seasonal and seasonal differences.
  if(tail(colnames(fit$var.coef),1)=="drift"){
    ARIMAModelInfo <- cbind(AR=fit$arma[1],Diff=fit$arma[6],MA=fit$arma[2],SeasonalAR=fit$arma[3],SeasonalDiff=fit$arma[7],SeasonalMA=fit$arma[4],Period=fit$arma[5],Drift=tail(fit$coef,1))
  }
  else {
    ARIMAModelInfo <- cbind(AR=fit$arma[1],Diff=fit$arma[6],MA=fit$arma[2],SeasonalAR=fit$arma[3],SeasonalDiff=fit$arma[7],SeasonalMA=fit$arma[4],Period=fit$arma[5],Drift=NA)
  }
  
  return (ARIMAModelInfo)
}
