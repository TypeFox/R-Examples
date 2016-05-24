#  LIMDIL.R

elda <- limdil <- function(response, dose, tested = rep(1, length(response)), group=rep(1,length(response)), observed = FALSE, confidence = 0.95, test.unit.slope = FALSE) 
#	Limiting dilution analysis
#	Gordon Smyth, Yifang Hu
#	21 June 2005. Last revised 18 August 2015.
{
	n <- length(response)
	if(n==0) stop("No data")
	if(length(dose) != n) stop("length(dose) doesn't match length(response)")
	if(length(tested) != n) {
		if(length(tested)==1)
			tested <- rep_len(tested,n)
		else
			stop("length(tested) doesn't match length(response)")
	}

#	Allow for structural zeros
	SZ <- response==0 & (dose==0 | tested==0)
	if(any(SZ)) {
		i <- !SZ
		out <- Recall(response=response[i],dose=dose[i],tested=tested[i],group=group[i],observed=observed,confidence=confidence,test.unit.slope=test.unit.slope)
		out$response <- response
		out$dose <- dose
		out$tested <- tested
		return(out)
	}

#	Check valid data
	y <- response/tested
	if (any(y < 0)) stop("Negative values for response or tested")
	if (any(y > 1)) stop("The response cannot be greater than the number tested")
	if (any(dose <= 0)) stop("dose must be positive")

	size <- 1 - confidence
	out <- list()
	f <- binomial(link = "cloglog")
	f$aic <- quasi()$aic

	group <- factor(group)
	num.group <- length(levels(group))
	groupLevel <- levels(group)

	out$response <- response
	out$tested <- tested
	out$dose <- dose
	out$group <- group
	out$num.group <- num.group
	class(out) <- "limdil"

	out$CI <- matrix(nrow=num.group,ncol=3)
	colnames(out$CI) <- c("Lower","Estimate","Upper")
	rownames(out$CI) <- paste("Group",levels(group))

#	Groupwise frequency estimates
	deviance0 <- dloglik.logdose <- FisherInfo.logdose <- dloglik.dose <- FisherInfo.dose <- 0
	for(i in 1:num.group) {
		index <- (group == groupLevel[i])
		fit0 <- eldaOneGroup(response=response[index],dose=dose[index],tested=tested[index],observed=observed,confidence=confidence,trace=FALSE)
		deviance0 <- deviance0 + fit0$deviance
		dloglik.logdose <- dloglik.logdose + fit0$dloglik.logdose
		FisherInfo.logdose <- FisherInfo.logdose + fit0$FisherInfo.logdose
		dloglik.dose <- dloglik.dose + fit0$dloglik.dose
		FisherInfo.dose <- FisherInfo.dose + fit0$FisherInfo.dose
		out$CI[i,] <- pmax(fit0$CI.frequency,1)
	}

#	Test for difference between groups
	if(num.group>1) {
		fitequal <- eldaOneGroup(response=response,dose=dose,tested=tested,observed=observed,confidence=confidence,trace=FALSE)
		dev.g <- pmax(fitequal$deviance - deviance0, 0)
		group.p <- pchisq(dev.g, df=num.group-1, lower.tail=FALSE)
		out$test.difference <- c(Chisq=dev.g, P.value=group.p, df=num.group-1)
	}

#	Test for unit slope
	if(test.unit.slope) {
		if(is.na(FisherInfo.logdose)) FisherInfo.logdose <- 0
		if(FisherInfo.logdose > 1e-15) {

#			Wald test
			if(num.group>1)
				fit.slope <- suppressWarnings(glm(y~group+log(dose), family=f, weights=tested))
			else
				fit.slope <- suppressWarnings(glm(y~log(dose), family=f, weights=tested))
			s.slope <- summary(fit.slope)
			est.slope <- s.slope$coef["log(dose)","Estimate"]
			se.slope <- s.slope$coef["log(dose)", "Std. Error"]
			z.wald <- (est.slope-1)/se.slope
			p.wald <- 2*pnorm(-abs(z.wald))
			out$test.slope.wald <- c("Estimate"=est.slope, "Std. Error"=se.slope, "z value"=z.wald, "Pr(>|z|)"=p.wald)

#			Likelihood ratio test
			dev <- pmax(deviance0 - fit.slope$deviance,0)
			z.lr <- sqrt(dev)*sign(z.wald)
			p.lr <- pchisq(dev, df = 1, lower.tail = FALSE)
			out$test.slope.lr <- c("Estimate"=NA, "Std. Error"=NA, "z value"=z.lr, "Pr(>|z|)"=p.lr)

#			Score tests for log(dose) and dose
			z.score.logdose <- dloglik.logdose / sqrt(FisherInfo.logdose)
			p.score.logdose <- 2*pnorm(-abs(z.score.logdose))
			z.score.dose <- dloglik.dose / sqrt(FisherInfo.dose)
			p.score.dose <- 2*pnorm(-abs(z.score.dose))
			out$test.slope.score.logdose <- c("Estimate"= NA, "Std. Error"=NA, "z value"=z.score.logdose,"Pr(>|z|)"=p.score.logdose)
			out$test.slope.score.dose <- c("Estimate"= NA, "Std. Error"=NA, "z value"=z.score.dose,"Pr(>|z|)"=p.score.dose)

		} else {
			out$test.slope.wald <- out$test.slope.lr <- out$test.slope.score.logdose <- out$test.slope.score.dose  <- c("Estimate"=NA, "Std. Error"=NA, "z value"=NA, "Pr(>|z|)"=1)
		}
	}

	out
}

print.limdil <- function(x, ...)
#	Print method for limdil objects
#	Yifang Hu and Gordon Smyth
#	20 February 2009. Last revised 31 January 2013.
{
	cat("Confidence intervals for frequency:\n\n")
	print(x$CI)
	
	if(!is.null(x$test.difference)) {
		difference <- x$test.difference
		cat("\nDifferences between groups:\n")
		cat("Chisq",difference[1], "on", difference[3], "DF, p-value:", format.pval(difference[2],4), "\n")
	}

	if(!is.null(x$test.slope.wald)) {
		a <- rbind(x$test.slope.wald, x$test.slope.lr, x$test.slope.score.logdose, x$test.slope.score.dose)
		a <- data.frame(a, check.names=FALSE)
		rownames(a) <- c("Wald test", "LR test", "Score test: log(Dose)", "Score test: Dose")
		cat("\nGoodness of fit (test log-Dose slope equals 1):\n")
		suppressWarnings(printCoefmat(a,tst.ind=1,has.Pvalue=TRUE,P.values=TRUE))
	}
}

plot.limdil <- function(x, col.group=NULL, cex=1, lwd=1, legend.pos="bottomleft", ...)
#	Plot method for limdil objects
#	Yifang Hu and Gordon Smyth
#	20 February 2009.  Last revised 6 February 2013.
{
	x$group <- factor(x$group)
	num.group <- nlevels(x$group)
	if(is.null(col.group)) 
		col.group <- 1:num.group
	else
		col.group <- rep(col.group,num.group)

	col <- x$group
	levels(col) <- col.group
	col <- as.character(col)
	dose <- x$dose
	maxx <- max(dose)	
	
	i <- x$response==x$tested
	x$response[i] <- x$response[i]-0.5

	nonres <- log(1-x$response/x$tested)
	if(num.group>1 && any(i)) nonres <- pmin(0,jitter(nonres))
		
	miny <- min(nonres)
	plot(x=1,y=1,xlim=c(0,maxx),ylim=c(min(miny,-0.5),0),xlab="dose (number of cells)",ylab="log fraction nonresponding",type="n",...)
	points(dose[!i],nonres[!i],pch=1,col=col[!i],cex=cex)
	points(dose[i],nonres[i],pch=6,col=col[i],cex=cex)

	for(g in 1:num.group) {
		abline(a=0,b=-1/x$CI[g,2],col=col.group[g],lty=1,lwd=lwd)
		abline(a=0,b=-1/x$CI[g,1],col=col.group[g],lty=2,lwd=lwd)
		abline(a=0,b=-1/x$CI[g,3],col=col.group[g],lty=2,lwd=lwd)
	}
	if(num.group>1) legend(legend.pos,legend=paste("Group",levels(x$group)),text.col=col.group,cex=0.6*cex)
	invisible(list(x=dose,y=nonres,group=x$group))
}

.limdil.allpos <- function(tested, dose, confidence, observed)
#	One-sided confidence interval when all assays are positive
#	Uses globally convergent Newton iteration
#	Yifang Hu.
#	Created 18 March 2009.  Last modified 18 Dec 2012.
{
	alpha <- 1 - confidence

	dosem <- min(dose)
	tested.group <- tested
	tested.sum <- sum(tested.group[dose == dosem])
	beta <- log(-log(1 - alpha^(1/tested.sum))) - log(dosem)

#	Starting value
	lambda <- exp(beta)
	if(observed) lambda <- -expm1(lambda)

#	Newton-iteration
	repeat {
		if(observed)
			f <- sum(tested*log(1-(1-lambda)^dose))-log(alpha)
		else
			f <- sum(tested*log(1-exp(-lambda*dose)))-log(alpha)
		if(observed)
			deriv <- sum(tested*(-dose)*(1-lambda)^(dose-1)/(1-(1-lambda)^dose)) 
		else
			deriv <- sum(tested*dose*exp(-dose*lambda)/(1-exp(-dose*lambda)))
		step <- f/deriv
		lambda <- lambda-step
		if(-step < 1e-6) break
	}
	lambda
}

eldaOneGroup <- function(response,dose,tested,observed=FALSE,confidence=0.95,tol=1e-8,maxit=100,trace=FALSE)
#	Estimate active cell frequency from LDA data
#	using globally convergent Newton iteration
#	Gordon Smyth
#	5 Dec 2012.  Last modified 30 Jan 2013.
{
	y <- response
	n <- tested
	d <- dose
	phat <- y/n
	size <- 1-confidence

#	Special case of all negative responses
	if(all(y < 1e-14)) {
	 	N <- sum(dose*tested)
		if (observed)
			U <- 1 - size^(1/N)
		else
			U <- -log(size)/N
		out <- list()
		out$CI.frequency <- c(Lower = Inf, Estimate = Inf, Upper = 1/U)
		out$deviance <- out$dloglik.logdose <- out$FisherInfo.logdose <- out$dloglik.dose <- out$FisherInfo.dose <- 0
		return(out)
	}

#	Special case of all positive responses
	if(all(phat > 1-1e-14)) {
		U <- .limdil.allpos(tested=tested,dose=dose,confidence=confidence,observed=observed)
		out <- list()
		out$CI.frequency <- c(Lower = 1/U, Estimate = 1, Upper = 1)
		out$deviance <- out$dloglik.logdose <- out$FisherInfo.logdose <- out$dloglik.dose <- out$FisherInfo.dose <- 0
		return(out)
	}

#	Starting value guaranteed to be left of the solution
	pmean <- mean(y)/mean(n)
	lambda <- -log1p(-pmean) / max(d)
	if(trace) cat(0,lambda,1/lambda,"\n")

#	Globally convergent Newton iteration
	iter <- 0
	repeat{
		iter <- iter+1
		if(iter > maxit) {
			warning("max iterations exceeded")
			break
		}
		p <- -expm1(-lambda*d)
		onemp <- exp(-lambda*d)

#		First derivative
		dloglik.lambda <- mean(n*d*(phat-p)/p)

#		Second derivative
		d2loglik.lambda <- -mean(n*phat*d*d*onemp/p/p)

#		Newton step
		step <- dloglik.lambda / d2loglik.lambda
		lambda <- lambda - step
		if(trace) cat(iter,lambda,1/lambda,step,"\n")
		if(abs(step) < tol) break
	}

#	Wald confidence interval for alpha
	alpha <- log(lambda)
	p <- -expm1(-lambda*d)
	onemp <- exp(-lambda*d)
	FisherInfo.alpha <- sum(n*d*d*onemp/p)*lambda^2
	SE.alpha <- 1/sqrt(FisherInfo.alpha)
	z <- qnorm( (1-confidence)/2, lower.tail=FALSE )
	CI.alpha <- c(Lower=alpha-z*SE.alpha,Estimate=alpha,Upper=alpha+z*SE.alpha)

#	Wald confidence interval for frequency
	if(observed)
		CI.frequency <- -1/expm1(-exp(CI.alpha))
	else
		CI.frequency <- exp(-CI.alpha)

#	Deviance
	f <- binomial(link="cloglog")
	deviance <- sum(f$dev.resid(phat,p,n))

#	Score test for log(dose) unit slope
	v <- p*onemp/n
	x <- log(d)
	eta <- alpha+x
	mu.eta <- f$mu.eta(eta)
	info.alpha <- mu.eta^2/v
	xmean <- sum(x*info.alpha)/sum(info.alpha)
	mu.beta <- (x-xmean)*mu.eta
	dloglik.beta <- sum(mu.beta*(phat-p)/v)
	FisherInfo.beta <- sum(mu.beta^2/v)
	z.scoretest <- dloglik.beta/sqrt(FisherInfo.beta)

#	Score test for dose
	x <- d
	xmean <- sum(x*info.alpha)/sum(info.alpha)
	mu.beta <- (x-xmean)*mu.eta
	dloglik.beta.dose <- sum(mu.beta*(phat-p)/v)
	FisherInfo.beta.dose <- sum(mu.beta^2/v)
	z.scoretest.dose <- dloglik.beta.dose/sqrt(FisherInfo.beta.dose)

	list(p=p,lambda=lambda,alpha=alpha,CI.alpha=CI.alpha,CI.frequency=CI.frequency,deviance=deviance,iter=iter,z.scoretest=z.scoretest,z.scoretest.dose=z.scoretest.dose,dloglik.logdose=dloglik.beta,FisherInfo.logdose=FisherInfo.beta,dloglik.dose=dloglik.beta.dose,FisherInfo.dose=FisherInfo.beta.dose)
}
