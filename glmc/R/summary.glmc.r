summary.glmc <- function(object, dispersion = NULL,
			correlation = FALSE, symbolic.cor = FALSE, ...)
{
    est.disp <- FALSE
    df.r <- object$df.residual

#    if(is.null(dispersion))	# calculate dispersion if needed
#        dispersion <- sum(object$weights*object$residuals^2)/object$df.r

        if(is.null(dispersion)){
	 dispersion <- {
	    if(any(object$family$family == c("poisson", "binomial")))  1
	    else if(df.r > 0) {
		est.disp <- TRUE
		if(any(object$weights==0))
		    warning("observations with zero weight ",
			    "not used for calculating dispersion")
		sum(object$weights*object$residuals^2)/ df.r
	    } else Inf
         }
      }

    ## calculate scaled and unscaled covariance matrix

    p <- object$rank
    if (p > 0) {
        p1 <- 1:p
        Qr <- object$qr
        aliased <- is.na(coef(object))  # used in print method
        ## WATCHIT! doesn't this rely on pivoting not permuting 1:p? -- that's quaranteed
        coef.p <- object$coefficients[Qr$pivot[p1]]
        if(!is.null(object$hessian)){
        covmat.unscaled<-object$hessian
        }else{
        covmat.unscaled <- chol2inv(Qr$qr[p1,p1,drop=FALSE])
        }

        dimnames(covmat.unscaled) <- list(names(coef.p),names(coef.p))
        covmat <- dispersion*covmat.unscaled
        var.cf <- diag(covmat)

        ## calculate coef table

        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err

        dn <- c("Estimate", "Std. Error")
        if(!est.disp) {
            pvalue <- 2*pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p),
                                         c(dn, "z value","Pr(>|z|)"))
        } else if(df.r > 0) {
            pvalue <- 2*pt(-abs(tvalue), df.r)
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p),
                                         c(dn, "t value","Pr(>|t|)"))
        } else { ## df.r == 0
            coef.table <- cbind(coef.p, Inf)
            dimnames(coef.table) <- list(names(coef.p), dn)
        }
        df.f <- NCOL(Qr$qr)
    } else {
        coef.table <- matrix(, 0, 4)
        dimnames(coef.table) <-
            list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0, 0)
        aliased <- is.na(coef(object))
        df.f <- length(aliased)
    }
    ## return answer

    ans <- c(object[c("call","terms","family","deviance", "aic",
		      "contrasts",
		      "df.residual","null.deviance","df.null","iter")],
	     list(deviance.resid = residuals(object, type = "deviance"),
		  coefficients = coef.table,
                  aliased = aliased,
		  dispersion = dispersion,
		  df = c(object$rank, df.r, df.f),
		  cov.unscaled = covmat.unscaled,
		  cov.scaled = covmat))

    if(correlation && p > 0) {
	dd <- sqrt(diag(covmat.unscaled))
	ans$correlation <-
	    covmat.unscaled/outer(dd,dd)
	ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.glmc"
    return(ans)
}
