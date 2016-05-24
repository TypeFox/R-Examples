FreqSurv <- function(Y, lin.pred, data)
{	

    ##
    y     <- as.vector(Y[,1])
    delta <- as.vector(Y[,2])
    Xmat  <- as.matrix(model.frame(lin.pred, data=data))
    ##
    fit.survreg <- survreg(as.formula(paste("Surv(y, delta) ", as.character(lin.pred)[1], as.character(lin.pred)[2])), dist="weibull", data=data)
    alpha    <- 1 / fit.survreg$scale
    
    ## log(kappa), log(alpha), log(beta)
    startVals <- c(-alpha*coef(fit.survreg)[1], log(alpha), -coef(fit.survreg)[-1]*alpha)
    ##
    fit0 <- nlm(logLike.weibull.Uni, p=startVals * runif(length(startVals), 0.9, 1.1),
                y=y, delta=delta, Xmat=Xmat,
                iterlim=1000, hessian=TRUE)
    ##
    if(fit0$code == 1 | fit0$code == 2)
    {
      myLabels <- c("log(kappa)", "log(alpha)", colnames(Xmat))
      value <- list(estimate=fit0$estimate, Finv=solve(fit0$hessian), logLike=-fit0$minimum, myLabels=myLabels)
      class(value) <- c("Freq", "Surv", "Ind", "WB")
      return(value)
    }


  ##
  invisible()
}