fp.order <- function(x, y, cox, gauss, xnames, ...)
{
#
# Provides ordering of input variables by LR test ranking via one step backward selection
#
    int <- as.numeric(!cox)
    nx <- ncol(x); nobs <- nrow(x)
    if(cox) {
        if(exists("coxph.fit")) fitter <- get("coxph.fit")
        else fitter <- getFromNamespace("coxph.fit","survival")
        fit <- fitter(x, y, ...)
        se <- sqrt(diag(fit$var))
		ns <- length(xnames)
        p.value <- numeric(ns)
        for(i in unique(xnames)) {
            ld <- sum(xnames==i)
            ll <- fitter(x[,xnames!=i, drop=FALSE], y, ...)$loglik[2]
            p.value[xnames==i] <- pchisq(2*(fit$loglik[2]-ll), df=ld, lower.tail=FALSE)
        }
        deviance <- -2 * fit$loglik
    }
    else {
        fit <- glm.fit(x, y, ...)    # full model
		df.r <- fit$df.residual
        dispersion <- if (gauss) 
						{
							if (df.r > 0) 
								sum(fit$residuals^2)/df.r
							else Inf
						}
						else 1
	    dev.full <- fit$deviance
	    dfs.full <- fit$df.residual
		ns <- length(xnames)
		p.value <- dev <- ld <- numeric(ns)
        for(i in unique(xnames)) {
            z <- glm.fit(x[,c(int,which(xnames!=i)+1), drop=FALSE], y, ...)
            dev[xnames==i] <- z$deviance
            ld[xnames==i] <- sum(xnames==i)
        }		
		dev <- if (gauss)  nobs * log(dev/nobs)  else  dev/dispersion
		dev.full <- if (gauss)  nobs * log(dev.full/nobs)  else  dev.full/dispersion
		dev <- pmax(0, dev - dev.full)
		p.value <- pchisq(dev, ld, lower.tail=FALSE)		
        deviance <- c(fit$null.deviance, fit$deviance)
	}
# dev.full linear
#
	x.order <- order(p.value)
#
    return(list(order = x.order, dev = deviance, df=c(nobs-int, nobs-nx)))
}
