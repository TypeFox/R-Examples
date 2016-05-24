summary.fit.linERR <-
function (object, ...) 
{
  call <- attr(object, "Call")
  coef.p <- attr(object,"details")$par
  var.cf <- diag(attr(object, "vcov"))
  dn <- c("Estimate", "Std. Error")
  s.err <- vector()
  tvalue <- vector()
  pvalue <- vector()
  
  names.covs <- c(colnames(attr(object, "covariates1")), attr(object, "covs_names"))
  if (!all(var.cf > 0))
  {
    warning("Hessian is not positive definite")
  }
  for (i in 1:length(coef.p))
  {
    if (var.cf[i] > 0)
    {
      s.err[i]  <- sqrt(var.cf[i])
      tvalue[i] <- coef.p[i]/s.err[i]
      pvalue[i] <- 2 * pnorm(-abs(tvalue[i]))
    }else{
      s.err[i]  <- NA
      tvalue[i] <- NA
      pvalue[i] <- NA  
    }
  }
    
  coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(names(coef.p), c(dn, "Test Stat.", "p-value"))
  rownames(coef.table) <- c("dose", names.covs)
  aic <- attr(object, "aic")
  dev <- attr(object, "deviance")
  inf.rsets <- length(attr(object, "rsets_2")) - 1
  ans <- list(coefficients=coef.table, aic=aic, dev=dev, call=call, inf.rsets=inf.rsets)
  
  class(ans) <- "summary.fit.linERR"
  return(ans)
}
