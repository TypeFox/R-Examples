"summary.frailtyPenal" <- function(object,level=.95, len=6, d=2, lab="hr", ...)
{

	x <- object
	if (!inherits(x, "frailtyPenal")) 
		stop("Object must be of class 'frailtyPenal'")
	
	nvar <- length(x$coef)
	
	if (is.null(x$coef)){
		cat("     Shared Gamma Frailty model: No covariates and no confidence interval\n")
	}else{
		if (nvar == 0){
			cat("No constant coefficients, only time-varying effects of the covariates \n")
		}else{
			z<-abs(qnorm((1-level)/2))
			co <- x$coef
			if(is.matrix(x$varH)){
				se <- sqrt(diag(x$varH))#[-1]
			}else{
				se <- sqrt(x$varH)
			}
			or <- exp(co)
			li <- exp(co-z * se)
			ls <- exp(co+z * se)
			r <- cbind(or, li, ls)
			dimnames(r) <- list(names(co), c(lab, paste(level*100,"%",sep=""), "C.I."))
			
			n<-r
		
			dd <- dim(n)
			n[n > 999.99] <- Inf
			a <- formatC(n, d, len,format="f")
	
			dim(a) <- dd
	
			if(length(dd) == 1){
				dd<-c(1,dd)
				dim(a)<-dd
				lab<-" "
			}else{
				lab <- dimnames(n)[[1]]
			}
	
			mx <- max(nchar(lab)) + 1
			cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
			
			for(i in (1):dd[1]) {
				lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
				cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
			}
		}
	}
}


