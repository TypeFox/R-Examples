##' @S3method cumres default
`cumres.default` <-
  function(model,score,information,
           residualfun,variable,
           data=model.frame(model),par=coef(model),
           R=1000, b=0, plots=min(R,50), seed=round(runif(1,1,1e9)),
           debug=FALSE, ...) {


    if(any(is.na(par))) stop("Over-parametrized model")
    ord <- order(variable)
    x <- variable[ord]
    n <- length(x)
    r <- residualfun(model)
    grad <- attributes(r)$grad    
    if (is.null(grad)) {
      if (!require("numDeriv")) stop("Supply gradient")
      grad <- numDeriv::jacobian(residualfun,par,...)
    }
    r <- r[ord]; grad <- grad[ord,,drop=FALSE]
    score <- score[ord,,drop=FALSE]
    Ii <-  solve(information)
    beta.iid <- Ii%*%t(score)
        
    hatW.MC <- function(x) {
      ord <- order(x)
      output <- .C("W",
                   R=as.integer(R), ## Number of realizations
                   b=as.double(b), ## Moving average parameter
                   n=as.integer(n), ## Number of observations
                   npar=as.integer(nrow(Ii)), ## Number of parameters (columns in design)
                   xdata=as.double(x), ## Observations to cummulate after 
                   rdata=as.double(r), ## Residuals (one-dimensional)
                   betaiiddata=as.double(beta.iid), ## Score-process
                   etarawdata=as.double(grad), ## Eta (derivative of terms in cummulated process W)
                   plotnum=as.integer(plots), ## Number of realizations to save (for later plot)
                   seed=as.integer(seed), ## Seed (will probably rely on R's rangen in future version)
                   KS=as.double(0), ## Return: Kolmogorov Smirnov statistics for each realization
                   CvM=as.double(0), ## Return: Cramer von Mises statistics for each realization
                   Wsd=as.double(numeric(n)), ## Return: Pointwise variance of W(x)
                   cvalues=as.double(numeric(R)), ## Return: value for each realization s.t.  +/- cvalue * Wsd contains W*
                   Ws=as.double(numeric(plots*n)), ## Return: Saved realizations (for plotting function)
                   Wobs=as.double(numeric(n)) ## Observed process
                   , PACKAGE="gof")
      return(list(output=output,x=x[ord]))
    }
        
    onesim <- hatW.MC(x)
    What <- matrix(onesim$output$Ws,nrow=n);
    ##    W <- cumsum(r[order(x0)]) 
    ##    matplot(x0,What,type="s", col="red", lty=1,pch=-1)
    ##    lines(onesim$output$Wobs ~ x0,type="s",lwd=2)

##    browser()
    
    res <- with(onesim$output,
                list(W=cbind(Wobs), What=list(What),
                     x=cbind(x),
                     KS=KS, CvM=CvM,
                     R=R, n=n, sd=cbind(Wsd), 
                     cvalues=cbind(cvalues), variable="x",
                     type="sem",
                     model=class(model)[1]) ##, onesim$output$WW)
                )
    class(res) <- "cumres"
    
    return(res)    
  }
