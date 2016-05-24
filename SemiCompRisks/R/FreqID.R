FreqID <- function(Y, lin.pred, data, model = "semi-Markov", frailty=TRUE)
{	

    ##
    y1     <- as.vector(Y[,1])
    delta1 <- as.vector(Y[,2])
    y2     <- as.vector(Y[,3])
    delta2 <- as.vector(Y[,4])
    Xmat1  <- as.matrix(model.frame(lin.pred[[1]], data=data))  
    Xmat2  <- as.matrix(model.frame(lin.pred[[2]], data=data))  
    Xmat3  <- as.matrix(model.frame(lin.pred[[3]], data=data))
    ##
    fit.survreg.1 <- survreg(as.formula(paste("Surv(y1, delta1) ", as.character(lin.pred[[1]])[1], as.character(lin.pred[[1]])[2])), dist="weibull", data=data)
    fit.survreg.2 <- survreg(as.formula(paste("Surv(y2, delta2) ", as.character(lin.pred[[2]])[1], as.character(lin.pred[[2]])[2])), dist="weibull", data=data)
    fit.survreg.3 <- survreg(as.formula(paste("Surv(y2, delta2) ", as.character(lin.pred[[3]])[1], as.character(lin.pred[[3]])[2])), dist="weibull", data=data)
    alpha1      <- 1 / fit.survreg.1$scale
    alpha2      <- 1 / fit.survreg.2$scale
    alpha3     	<- 1 / fit.survreg.3$scale
    
    
    startVals     <- c(-alpha1*coef(fit.survreg.1)[1], log(alpha1),
                       -alpha2*coef(fit.survreg.2)[1], log(alpha2),
                       -alpha3*coef(fit.survreg.3)[1], log(alpha3))
    if(frailty == TRUE) startVals <- c(startVals, 0.5)
    startVals     <- c(startVals,
                       -coef(fit.survreg.1)[-1] * alpha1,
                       -coef(fit.survreg.2)[-1] * alpha2,
                       -coef(fit.survreg.3)[-1] * alpha3)
    
    if(model == "Markov")
    {
        ##
        fit0  <- nlm(logLike.weibull.SCR, p=startVals * runif(length(startVals), 0.9, 1.1),
        y1=y1, delta1=delta1, y2=y2, delta2=delta2, Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,
        iterlim=1000, hessian=TRUE)
        
        fit1  <- nlm(logLike.weibull.SCR, p=runif(length(startVals), -1, 1),
        y1=y1, delta1=delta1, y2=y2, delta2=delta2, Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,
        iterlim=1000, hessian=TRUE)
    }
    if(model == "semi-Markov")
    {
        ##
        fit0  <- nlm(logLike.weibull.SCR.SM, p=startVals * runif(length(startVals), 0.9, 1.1),
        y1=y1, delta1=delta1, y2=y2, delta2=delta2, Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,
        iterlim=1000, hessian=TRUE)
        
        fit1  <- nlm(logLike.weibull.SCR.SM, p=runif(length(startVals), -1, 1),
        y1=y1, delta1=delta1, y2=y2, delta2=delta2, Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, frailty=frailty,
        iterlim=1000, hessian=TRUE)
    }

    ##
    if(fit0$code == 1 | fit0$code == 2)
    {
      ##
      myLabels <- c("log(kappa1)", "log(alpha1)", "log(kappa2)", "log(alpha2)", "log(kappa3)", "log(alpha3)")
      if(frailty == TRUE) myLabels <- c(myLabels, "log(theta)")
      myLabels <- c(myLabels, colnames(Xmat1), colnames(Xmat2), colnames(Xmat3))
      nP <- c(ncol(Xmat1), ncol(Xmat2), ncol(Xmat3))
      ##
      value <- list(estimate=fit0$estimate, Finv=solve(fit0$hessian), logLike=-fit0$minimum, myLabels=myLabels, frailty=frailty, nP=nP, Xmat=list(Xmat1, Xmat2, Xmat3))
      
      if(model == "semi-Markov")
      {
          class(value) <- c("Freq", "ID", "Ind", "WB", "semi-Markov")
      }
      if(model == "Markov")
      {
          class(value) <- c("Freq", "ID", "Ind", "WB", "Markov")
      }
      
      return(value)
    }
    
  ##
  invisible()
}
