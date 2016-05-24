"nlsBoot"<-function(nls, niter=999){

	if (!inherits(nls, "nls"))
		stop("Use only with 'nls' objects")

	data2 <- eval(nls$data, sys.frame(0))
	fitted1 <- fitted(nls)
	resid1 <- resid(nls)
	var1 <- all.vars(formula(nls)[[2]])
	
	l1 <- lapply(1:niter, function(i){
		data2[,var1] <- fitted1 + sample(scale(resid1, scale=FALSE), replace=TRUE);
		nls2 <- try(update(nls, start=as.list(coef(nls)), data=data2), silent=TRUE);
		if(inherits(nls2, "nls"))
			return(list(coef=coef(nls2), rse=summary(nls2)$sigma))
		})

	if(sum(sapply(l1, is.null)) > niter/2) stop(paste("Procedure aborted: the fit only converged in", round(sum(sapply(l1, is.null))/niter), "% during bootstrapping"))

	tabboot <- sapply(l1[!sapply(l1, is.null)], function(z) z$coef)
	rseboot <- sapply(l1[!sapply(l1, is.null)], function(z) z$rse)
	recapboot <- t(apply(tabboot, 1, quantile, c(.5, .025, .975))); colnames(recapboot) <- c("Median", "2.5%", "97.5%")
	estiboot <- t(apply(tabboot, 1, function(z) c(mean(z), sd(z)))); colnames(estiboot) <- c("Estimate", "Std. error")
	
	serr <- sum(sapply(l1, is.null))
	if(serr > 0) warning(paste("The fit did not converge", serr, "times during bootstrapping"))
	
	listboot <- list(coefboot = t(tabboot), rse = rseboot, bootCI = recapboot, estiboot = estiboot)
	class(listboot) <- "nlsBoot"
	return(listboot)
	
}


"plot.nlsBoot"<-function(x, type=c("pairs","boxplot"), mfr=c(ceiling(sqrt(ncol(x$coefboot))),ceiling(sqrt(ncol(x$coefboot)))),ask=FALSE, ...){
	if (!inherits(x, "nlsBoot"))
		stop("Use only with 'nlsBoot' objects")
	tab <- x$coefboot
	np <- ncol(tab)
 	def.par <- par(no.readonly = TRUE)	
	if(type[1] == "pairs"){
		if(ask) par(ask=TRUE, mar=c(4,4,3,1))
		if(!ask){
			lay <- lower.tri(matrix(0,(np-1),(np-1)), TRUE)
			lay[which(lay, TRUE)] <- 1:choose(np,2)
			layout(lay)
			par(mar=c(5,4,0.2,0.2))
		}
		for(i in 1:(np-1))
			for(j in (i+1):np)
				plot(tab[,i], tab[,j], xlab=colnames(tab)[i], ylab=colnames(tab)[j], pch="+")
	}
	if(type[1] == "boxplot"){ 
		if(ask) par(ask=TRUE, mar=c(4,4,3,1))
		if(!ask) par(mfrow=mfr, mar=c(4,4,3,1))
		for(i in 1:np){
			boxplot(tab[,i],main=colnames(tab)[i])
		}
	}
	par(def.par)
}


"print.nlsBoot" <- function (x, ...) {
	if (!inherits(x, "nlsBoot"))
		stop("Use only with 'nlsBoot' objects")
	cat("Bootstrap resampling\n")
	cat("\n")
	sumry <- array("", c(1, 4), list(1:1, c("vector", "length", "mode", "content")))
	sumry[1, ] <- c("$rse", length(x$rse), mode(x$rse), "Bootstrap residual errors")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
	sumry <- array("", c(3, 4), list(1:3, c("data.frame", "nrow", "ncol", "content")))
	sumry[1, ] <- c("$coefboot", nrow(x$coefboot), ncol(x$coefboot), "Bootstrap parameter estimates")
	sumry[2, ] <- c("$estiboot", nrow(x$estiboot), ncol(x$estiboot), "Bootstrap estimates and std. error")
	sumry[3, ] <- c("$bootCI", nrow(x$bootCI), ncol(x$bootCI), "Bootstrap medians and 95% CI")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}

"summary.nlsBoot" <- function (object, ...) {
	if (!inherits(object, "nlsBoot"))
		stop("Use only with 'nlsBoot' objects")
	cat("\n------\n")
	cat("Bootstrap statistics\n")
	print(object$estiboot)
	cat("\n------\n")
	cat("Median of bootstrap estimates and percentile confidence intervals\n")
	print(object$bootCI)
	cat("\n")
}

