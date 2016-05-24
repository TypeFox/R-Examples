"nlsJack"<-function(nls){
	if (!inherits(nls, "nls"))
		stop("Use only with 'nls' objects")

	c1 <- nls$call
	data <- eval(c1$data, sys.frame(0))
	c1$start <- as.list(coef(nls))
	nl <- nrow(data)
	np <- length(coef(nls))
	
	l1 <- lapply(1:nl, function(z){
		nls2 <- update(nls, data=data[-z,], start=as.list(coef(nls)));
		return(list(coef=coef(nls2), sigma=summary(nls2)$sigma, rss=sum(residuals(nls2)^2), dfb=abs(coef(nls2)-coef(nls))/(summary(nls)$parameters[,2])))
	})

	tabjack <- t(sapply(l1, function(z) z$coef))
	rsejack <- sapply(l1, function(z) z$sigma)
	rssjack <- sapply(l1, function(z) z$rss)
	dfb <- t(sapply(l1, function(z) z$dfb))
	reldif	<- apply(tabjack, 1, function(x) 100*abs(x-coef(nls))/coef(nls))
	pseudo <- t(apply((nl-1) * tabjack, 1, function(z) nl * coef(nls) - z))
	estijack <- cbind.data.frame(Estimates=colSums(pseudo) / nl, Bias=coef(nls) - colSums(pseudo) / nl)
	sum1 <- crossprod(t(t(pseudo) - estijack$Estimates))
	varjack <- (1 / (nl * (nl - 1))) * sum1
	student95 <- qt(0.975, df = nl - np)
	ICjack <- cbind.data.frame(Esti = estijack$Estimates, Low = estijack$Estimates - student95 * sqrt(diag(varjack)), Up = estijack$Estimates + student95 * sqrt(diag(varjack)))

	listjack	<-list(estijack=estijack, coefjack=tabjack, reldif=reldif, rse=rsejack, rss=rssjack ,dfb=dfb, jackCI=ICjack)
	class(listjack) <- "nlsJack"
	return(listjack)
}

"plot.nlsJack"<-function(x, mfr=c(nrow(x$reldif),1), ask=FALSE, ...){
	if (!inherits(x, "nlsJack"))
		stop("Use only with 'nlsJack' objects")
	if(ask) par(ask=TRUE,mar=c(4,4,3,1))
	if(!ask) par(mfrow=mfr,mar=c(4,4,3,1))
	for(i in 1:nrow(x$reldif)){
		plot(x$reldif[i,],type="h",main=rownames(x$reldif)[i],xlab="Observation #",ylab="Rel Diff (%)",ylim=c(0,1.2*max(x$reldif[i,])))
		for(j in 1:nrow(x$dfb)){
			if(x$dfb[j,i]>(2/sqrt(nrow(x$dfb))))
				text(j,x$reldif[i,j],"*",col="red",cex=2)
		}
	}
	par(mfrow=c(1,1),ask=FALSE)
}

"print.nlsJack" <- function (x, ...) {
	if (!inherits(x, "nlsJack"))
		stop("Use only with 'nlsJack' objects")
	cat("Jackknife resampling\n")
	cat("\n")
	sumry <- array("", c(3, 4), list(1:3, c("vector", "length", "mode", "content")))
	sumry[1, ] <- c("$estijack", length(x$estijack), mode(x$estijack), "jackknife estimates and bias")
	sumry[2, ] <- c("$rse", length(x$rse), mode(x$rse), "residual errors")
	sumry[3, ] <- c("$rss", length(x$rss), mode(x$rss), "residual sum of squares")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
	sumry <- array("", c(4, 4), list(1:4, c("data.frame", "nrow", "ncol", "content")))
	sumry[1, ] <- c("$coefjack", nrow(x$coefjack), ncol(x$coefjack), "parameter estimates")
	sumry[2, ] <- c("$reldif", nrow(x$reldif), ncol(x$reldif), "relative diff. of parameter estimates")
	sumry[3, ] <- c("$dfb", nrow(x$dfb), ncol(x$dfb), "DF-beta")
	sumry[4, ] <- c("$jackCI", nrow(x$jackCI), ncol(x$jackCI), "jackknife confidence intervals")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}

"summary.nlsJack" <- function (object, ...) {
	if (!inherits(object, "nlsJack"))
		stop("Use only with 'nlsJack' objects")
	cat("\n------\n")
	cat("Jackknife statistics\n")
	print(object$estijack)
	cat("\n------\n")
	cat("Jackknife confidence intervals\n")
	print(object$jackCI[,c("Low","Up")])
	cat("\n------\n")
	cat("Influential values\n")
	inf1 <- which(object$dfb>(2/sqrt(nrow(object$dfb))), arr.ind=TRUE)
	for(i in 1:nrow(inf1))
		cat("* Observation", inf1[i,1], "is influential on", rownames(object$estijack)[inf1[i,2]], "\n")
	cat("\n")
}
