"nlsResiduals" <- function(nls){
	if (!inherits(nls, "nls"))
		stop("Use only with 'nls' objects")

	#which = 1 residuals
	resi1		<- cbind(fitted(nls), residuals(nls))
	colnames(resi1)	<- c("Fitted values", "Residuals")

	#which = 2 residuals
	ecartres	<- summary(nls)$sigma
	nresiduals 	<- (residuals(nls)-mean(residuals(nls)))/ecartres
	std95		<- qt(0.975,df=length(resid(nls))-length(coef(nls)))
	resi2		<- cbind(fitted(nls), nresiduals)
	colnames(resi2)	<- c("Fitted values", "Standardized residuals")

	#which = 3 residuals
	ecartres	<- summary(nls)$sigma
	nresiduals	<- residuals(nls)/ecartres
	resi3		<- cbind(fitted(nls), sqrt(abs(nresiduals)))
	colnames(resi3)	<- c("Fitted values", "Sqrt abs. standardized residuals")

	#which = 4 residuals
	resiminus	<-vector()
	resiplus	<-vector()
       	for(i in 1:(length(residuals(nls))-1)){
               	resiminus[i]	<- residuals(nls)[i]
               	resiplus[i]	<- residuals(nls)[i+1]
        }
	resi4	<- cbind(resiminus, resiplus)
	colnames(resi4)	<- c("Residuals i", "Residuals i+1")

	listresi	<- list(resi1=resi1, resi2=resi2, resi3=resi3, resi4=resi4, std95=std95)	
	class(listresi)	<- "nlsResiduals"
	return(listresi)
}

"plot.nlsResiduals" <- function(x, which=0, ...){

	"hist.nlsResiduals" <- function(x, ...){
	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")
	hist(x$resi1[,"Residuals"], main="Residuals", xlab="Residuals")
	}

	"boxplot.nlsResiduals" <- function(x, ...){
	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")
	boxplot(x$resi1[,"Residuals"], main="Residuals")
	}

	"qq.nlsResiduals" <- function(x){
	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")

	qqploty	<- sort(x$resi2[,2])
	qqplotx	<- qnorm(ppoints(nrow(x$resi2)))
	qy	<- quantile(qqploty,c(0.25,0.75))
	qx	<- qnorm(c(0.25,0.75))
	slope	<- diff(qy)/diff(qx)
	ori	<- qy[1]-slope*qx[1]     
	
	qqnorm(x$resi2[,2], main="Normal Q-Q Plot of\n Standardized Residuals")
	qqline(x$resi2[,2])
	}

	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")
	if(!(which %in% 0:6)) 
		stop("\n Expected 'which':\n 1 = un-transformed\n 2 = normed\n 3 = sqrt absolute normed\n 4 = auto-correlation\n 5 = histogram\n 6 = qq-plot")
	if(which == 0){
		def.par <- par(no.readonly = TRUE)
        par(mfrow=c(2,2))
        plot.nlsResiduals(x, which=1)
        plot.nlsResiduals(x, which=2)
        plot.nlsResiduals(x, which=4)
        qq.nlsResiduals(x)
        par(def.par)
	}
	if(which == 1){
		plot(x$resi1, xlab="Fitted values", ylab="Residuals", main="Residuals")
		abline(h=0, lty=2)
	}
	if(which == 2){
		yrange	<- c(min(min(x$resi2[,2]), -x$std95), max(max(x$resi2[,2]), x$std95))
		plot(x$resi2, xlab="Fitted values", ylab="Standardized residuals", ylim=yrange, main="Standardized Residuals")
		abline(h=0,lty=2); abline(h=x$std95); abline(h=-x$std95)
	}
	if(which == 3){
		plot(x$resi3, xlab="Fitted values", ylab=expression(sqrt(abs("Standardized residuals"))), main="Sqrt abs residuals")
	}
	if(which == 4){
		plot(x$resi4, xlab="Residuals i", ylab="Residuals i+1", main="Autocorrelation")
		abline(h=0, lty=2)
	}
	if(which == 5){
		hist(x)
	}
	if(which == 6){
		qq.nlsResiduals(x)
	}
}


"test.nlsResiduals" <- function(x){
	"runs.test" <- function (x, alternative = c("two.sided", "less", "greater")){
		if(!is.factor(x))
			stop("x is not a factor")
		if(any(is.na(x)))
			stop("NAs in x")
		if(length(levels(x)) != 2)
			stop("x does not contain dichotomous data")
		alternative <- match.arg(alternative)
		DNAME <- deparse(substitute(x))
		n <- length(x)
		R <- 1 + sum(as.numeric(x[-1] != x[-n]))
		n1 <- sum(levels(x)[1] == x)
		n2 <- sum(levels(x)[2] == x)
		m <- 1 + 2*n1*n2 / (n1+n2)
		s <- sqrt(2*n1*n2 * (2*n1*n2 - n1 - n2) / ((n1+n2)^2 * (n1+n2-1)))
		STATISTIC <- (R - m) / s
		METHOD <- "Runs Test"
		if(alternative == "two.sided")
			PVAL <- 2 * pnorm(-abs(STATISTIC))
		else if(alternative == "less")
			PVAL <- pnorm(STATISTIC)
		else if(alternative == "greater")
			PVAL <- pnorm(STATISTIC, lower.tail = FALSE)
		else stop("irregular alternative")
		names(STATISTIC) <- "Standard Normal"
		structure(list(statistic = STATISTIC,
			alternative = alternative,
			p.value = PVAL,
			method = METHOD,
			data.name = DNAME),
			class = "htest")
	}

	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")
	#Shapiro-Wilk test
	stdres <- x$resi2[,2]
	shapi	<- shapiro.test(stdres)
	cat("\n------")
	print(shapi)

	#Test runs
	run	<- vector(length=nrow(x$resi1))
	for(i in 1:nrow(x$resi1)){
		if(x$resi2[i,2]<0) run[i] <- "N" 
        if(x$resi2[i,2]>0) run[i] <- "P"
	}
	cat("\n------")
	runs.test(as.factor(run))
}

"print.nlsResiduals" <- function(x, ...){
	if (!inherits(x, "nlsResiduals"))
		stop("Use only with 'nlsResiduals' objects")
	cat("Residuals\n")
	cat("\n")
	sumry <- array("", c(1, 4), list(1:1, c("vector", "length", "mode", "content")))
	sumry[1, ] <- c("$std95", length(x$std95), mode(x$std95), "Student value for alpha = 0.05")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
	sumry <- array("", c(4, 4), list(1:4, c("matrix", "nrow", "ncol", "content")))
	sumry[1, ] <- c("$resi1", nrow(x$resi1), ncol(x$resi1), "fitted values vs. non-transf. resid")
	sumry[2, ] <- c("$resi2", nrow(x$resi2), ncol(x$resi2), "fitted values vs. standardized resid")
	sumry[3, ] <- c("$resi3", nrow(x$resi3), ncol(x$resi3), "fitted values vs. sqrt abs std resid")
	sumry[4, ] <- c("$resi4", nrow(x$resi4), ncol(x$resi4), "resid i vs. resid i+1")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}
