"preview" <- function(formula, data, start, variable=1){

	"formula2function"<-function(formu){
		arg1		<- all.vars(formu)
		arg2		<- vector("list",length(arg1))
		names(arg2)	<- arg1
		Args		<- do.call("alist",arg2)
		fmodele		<- as.function(c(Args,formu))
		return(fmodele)
	}

	f1 <- formula2function(formula(formula)[[3]])
	vardep <- all.vars(formula[[2]])
	varindep <- intersect(all.vars(formula[[3]]), colnames(data))
	predic <- do.call(f1, as.list(c(start, data[varindep])))
	rss1 <- signif(sum((predic-data[vardep])^2), 3)

	plot(data[c(variable, which(colnames(data)==vardep))], ylab="Predicted", main="", ylim=c(min(data[vardep],predic), max(data[vardep],predic)))
	points(cbind.data.frame(data[variable], predic), pch="+", col="red")
	cat("\nRSS: ", rss1,"\n")

}

"plotfit" <- function(x, smooth=FALSE, variable=1, xlab=NULL, ylab=NULL, pch.obs=1, pch.fit="+", lty=1, lwd=1, col.obs="black", col.fit="red", ...){
	if (!inherits(x, "nls"))
		stop("Use only with 'nls' objects")
	d <- eval(x$call$data, sys.frame(0))
	vardep <- all.vars(formula(x)[[2]])
	varindep <- intersect(all.vars(formula(x)[[3]]), colnames(d))
	variable1 <- which(varindep == colnames(d)[variable])
	if (smooth & length(varindep)!=1) 
        	stop("smooth option is only possible when the number of independent variables equals 1")
	if(smooth | smooth=="T"){
        w0 <- list(seq(min(d[,varindep]), max(d[,varindep]), len=1000))
		names(w0) <- varindep
		if(is.null(xlab)) xlab <- varindep
		if(is.null(ylab)) ylab <- vardep
		plot(d[c(varindep, vardep)], xlab=xlab, ylab=ylab, pch=pch.obs, col=col.obs, ...)
		lines(w0[[1]], predict(x,new=w0), col=col.fit, lty=lty, lwd=lwd)
	}
	else{
		if(is.null(xlab)) xlab <- varindep[variable1]
		if(is.null(ylab)) ylab <- vardep
		plot(d[,vardep] ~ d[,varindep[variable1]], xlab=xlab, ylab=ylab, pch=pch.obs, col=col.obs, ...)
		points(d[,varindep[variable1]], predict(x), pch=pch.fit, col=col.fit)
	}
}

"overview" <- function(x){
	if (!inherits(x, "nls"))
		stop("Use only with 'nls' objects")
	cat("\n------")
	print(summary(x))
	cat("------\n")
	cat("Residual sum of squares:", signif(sum(residuals(x)^2), 3),"\n\n")
	n <- length(residuals(x))
	np <- length(coef(x))
	esti <- summary(x)$parameters[,"Estimate"]
	ster <- summary(x)$parameters[,"Std. Error"]
	t95 <- qt(0.975, df=(n-np))
	binf <- esti - t95 * ster
	bsup <- esti + t95 * ster
	cat("------\n")
	cat("t-based confidence interval:\n")
	print(cbind.data.frame("2.5%" = binf, "97.5%" = bsup))
	cat("\n")
	cat("------\n")
	cat("Correlation matrix:\n")
	print(summary(x, correlation = TRUE)$correlation)
	cat("\n")
}
