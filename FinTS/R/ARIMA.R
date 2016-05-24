ARIMA <- function(x, order = c(0, 0, 0), 
      seasonal = list(order = c(0, 0, 0), period = NA),
      xreg = NULL, include.mean = TRUE, transform.pars = TRUE,
      fixed = NULL, init = NULL, method = c("CSS-ML", "ML", "CSS"),
      n.cond, optim.control = list(), kappa = 1e6, Box.test.lag=NULL,
      Box.test.df = c("net.lag", "lag"), 
      type = c("Ljung-Box", "Box-Pierce", "rank")){
##
## 1.  arima
##
  fit <- arima(x=x, order=order, seasonal=seasonal, xreg=xreg,
               include.mean=include.mean, transform.pars=transform.pars, 
               fixed=fixed, init=init, method=method, n.cond=n.cond,
               optim.control=optim.control, kappa=kappa)
##
## 2.  Compute desired number of lags and degrees of freedom
##     for Box.test 
##
#  2.1.  number of parameters estimated (apart from 'intercept')
  vc <- vcov(fit)
# CAN NOT use coef(fit) here,
# because it includes parameters fixed as well as estimated.
# vcov includes only parameters ESTIMATED  
  int <- ('intercept' %in% dimnames(vc)[[1]])
  kPars <- (dim(vc)[1] - int)
# 2.2.  Box.test.lag 
  if(is.null(Box.test.lag))
    Box.test.lag <- round(log(sum(!is.na(x))))  
  Lag <- max(kPars+1, Box.test.lag)
# 2.3.  Box.test.df
  df. <- Box.test.df
  {
    if(is.character(df.)){
      Box.test.df <- match.arg(Box.test.df)
      df. <- (Lag - (Box.test.df == "net.lag") * kPars)
    }
    else
      if(df.<=0) df. <- 1
  }
# 2.4.  'Ljung-Box' or 'Box-Pierce'?    
  tp <- match.arg(type)
##
## 3.  Compute AutocorTest 
##  
  LjB <- AutocorTest(fit$resid, lag=Lag, type=tp, df=df.)
  fit$Box.test <- LjB
##
## 4.  'xreg'?  
##
  if(!is.null(xreg)){
    varX <- var(xreg)
    k <- dim(varX)[1]
    if(length(k)==0)k <- 1 
    b <- coef(fit)
    bx <- b[length(b)-(k-1):0]
    var.expl <- crossprod(varX %*% bx, bx)
    fit$r.squared <- min(1, var.expl/var(x))
  }
##
## 5.  Done
##  
  fit
}
