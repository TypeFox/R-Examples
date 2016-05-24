"summary.ahaz" <- function(object, ...)
  {
    ## Purpose: Summary of 'ahaz' object
    ## ----------------------------------------------------------------------
    ## Arguments:
    ##   object     : 'ahaz' object
    ##     ...      : additional arguments
    ## ----------------------------------------------------------------------
    ## Author: Anders Gorst-Rasmussen
    
   
    if (!object$univar)
      {
         # Try to handle situation w/collinearity
         # OK if no collinearity, somewhat slow otherwise
        Dinv<-matrix(NA, nrow = object$nvars, ncol = object$nvars)
        beta<-rep(NA, object$nvars)
        
        Dqr <- qr(object$D)
        use<-Dqr$pivot[1:Dqr$rank]
        if(Dqr$rank == object$nvars)
          Dinv<-qr.solve(Dqr)
        else
          Dinv[use,use]<-solve(object$D[use, use])
        beta[use]<-Dinv[use, use] %*% object$d[use]
                                            
        cov <- (Dinv %*% object$B %*% Dinv) 
        se <- diag(cov) ^ .5
        if(Dqr$rank==object$nvars)
          wtest<-drop(object$d %*% solve(object$B) %*%object$d )
        else
          wtest<-NA
        waldtest<-c("test"=wtest, "pvalue"=1 - pchisq(wtest,df = object$nvars),"df"=object$nvars)
      }
    else{
      beta<- object$d/object$D
      cov <- object$B / object$D^2 
      se <- cov^.5
      waldtest<-NULL
    }
    z <- beta / se
    pval <- 2 * (1 - ifelse(z > 0, pnorm(z), pnorm(-z)))
    coeff <- cbind(beta, se, z, pval)

    colnames(coeff) <- c("Estimate",
                         ifelse(object$robust,"Robust SE","Std. Error"),
                         "Z value", "Pr(>|z|)")

    if(is.null(object$data$colnames))
      rownames(coeff) <- paste("Coeff.",1:length(beta))
    else
      rownames(coeff) <- object$data$colnames
    
    out <- list("call" = object$call,"coefficients" = coeff,"cov"=cov, "nobs" = object$nvars,
                "nvars" = object$nobs,"waldtest"=waldtest,"univar" = object$univar)
    class(out) <- "summary.ahaz"
    return(out)
  }
