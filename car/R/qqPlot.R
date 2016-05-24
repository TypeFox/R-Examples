# Quantile-comparison plots (J. Fox)

# last modified 30 September 2009 by J. Fox
# November 2009 by S. Weisberg -- changed to use showLabels for point identification
# 14 April 2010: set id.n = 0. J. Fox
# 1 June 2010: set reps=100 in qqPlot.lm. J. Fox
# 28 June 2010: fixed labeling bug S. Weisberg
# 11 March 2011: moved up ... argument. J. Fox
# 23 May 2012: line="none" now honored in qqPlot.default. J. Fox
# 2 May 2013: qqPlot.lm() now works with "aov" objects (fixing problem reported by Thomas Burk). J. Fox
# 2015-12-12: allow vectorized col, pch, and cex arguments (suggestion of Emmanuel Curis)

qqp <- function(...) qqPlot(...)

qqPlot<-function(x, ...) {
	UseMethod("qqPlot")
}

qqPlot.default <- function(x, distribution="norm", ..., ylab=deparse(substitute(x)),
		xlab=paste(distribution, "quantiles"), main=NULL, las=par("las"),
		envelope=.95,  
		col=palette()[1], col.lines=palette()[2], lwd=2, pch=1, cex=par("cex"), 
		line=c("quartiles", "robust", "none"), 
		labels = if(!is.null(names(x))) names(x) else seq(along=x),
		id.method = "y", 
		id.n = if(id.method[1]=="identify") Inf else 0,
		id.cex=1, id.col=palette()[1], grid=TRUE)
{
	line <- match.arg(line)
	good <- !is.na(x)
	ord <- order(x[good])
	if (length(col) == length(x)) col <- col[good][ord]
	if (length(pch) == length(x)) pch <- pch[good][ord]
	if (length(cex) == length(x)) cex <- cex[good][ord]
	ord.x <- x[good][ord]
	ord.lab <- labels[good][ord]
	q.function <- eval(parse(text=paste("q", distribution, sep="")))
	d.function <- eval(parse(text=paste("d", distribution, sep="")))
	n <- length(ord.x)
	P <- ppoints(n)
	z <- q.function(P, ...)
	plot(z, ord.x, type="n", xlab=xlab, ylab=ylab, main=main, las=las)
	if(grid){
		grid(lty=1, equilogs=FALSE)
		box()}
	points(z, ord.x, col=col, pch=pch, cex=cex)
	if (line == "quartiles" || line == "none"){
		Q.x <- quantile(ord.x, c(.25,.75))
		Q.z <- q.function(c(.25,.75), ...)
		b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
		a <- Q.x[1] - b*Q.z[1]
		if (line == "quartiles") abline(a, b, col=col.lines, lwd=lwd)
	}
	if (line=="robust") {
		coef <- coef(rlm(ord.x ~ z))
		a <- coef[1]
		b <- coef[2]
		abline(a, b, col=col.lines, lwd=lwd)
	}
	conf <- if (envelope == FALSE) .95 else envelope
	zz <- qnorm(1 - (1 - conf)/2)
	SE <- (b/d.function(z, ...))*sqrt(P*(1 - P)/n)
	fit.value <- a + b*z
	upper <- fit.value + zz*SE
	lower <- fit.value - zz*SE
	if (envelope != FALSE) {
		lines(z, upper, lty=2, lwd=lwd, col=col.lines)
		lines(z, lower, lty=2, lwd=lwd, col=col.lines)
	}
	showLabels(z, ord.x, labels=ord.lab,
			id.method = id.method, id.n = id.n, id.cex=id.cex, id.col=id.col)
}

qqPlot.lm <- function(x, xlab=paste(distribution, "Quantiles"),
		ylab=paste("Studentized Residuals(", deparse(substitute(x)), ")", sep=""), main=NULL,
		distribution=c("t", "norm"), line=c("robust", "quartiles", "none"), las=par("las"),
		simulate=TRUE, envelope=.95,  reps=100, 
		col=palette()[1], col.lines=palette()[2], lwd=2, pch=1, cex=par("cex"),
		labels, id.method = "y", 
		id.n = if(id.method[1]=="identify") Inf else 0, id.cex=1, 
		id.col=palette()[1], grid=TRUE, ...){
	result <- NULL
	distribution <- match.arg(distribution)
	line <- match.arg(line)
	rstudent <- rstudent(x)
	if (missing(labels)) labels <- names(rstudent)
	good <- !is.na(rstudent)
	rstudent <- rstudent[good]
	labels <- labels[good]
	sumry <- summary.lm(x)
	res.df <- sumry$df[2]
	if(!simulate)
		result <- qqPlot(rstudent, distribution=if (distribution == "t") "t" else "norm", df=res.df-1, line=line,
				main=main, xlab=xlab, ylab=ylab, las=las, envelope=envelope, 
				col=col, col.lines=col.lines, lwd=lwd, pch=pch, cex=cex,
				labels=labels, id.method=id.method, id.n=id.n, id.cex=id.cex,
				id.col=id.col, ...)
	else {
		n <- length(rstudent)        
		ord <- order(rstudent)
		ord.x <- rstudent[ord]
		ord.lab <- labels[ord]
		P <- ppoints(n)
		z <- if (distribution == 't') qt(P, df=res.df-1) else qnorm(P)
		plot(z, ord.x, type="n", xlab=xlab, ylab=ylab, main=main, las=las)
		if(grid) grid(lty=1, equilogs=FALSE)
		points(z, ord.x, pch=pch, col=col, cex=cex)
		yhat <- na.omit(fitted.values(x))
		S <- sumry$sigma
		Y <- matrix(yhat, n, reps) + matrix(rnorm(n*reps, sd=S), n, reps)
		X <- model.matrix(x)
		rstud <- apply(rstudent(lm(Y ~ X - 1)), 2, sort)
		lower <- apply(rstud, 1, quantile, prob=(1 - envelope)/2)
		upper <- apply(rstud, 1, quantile, prob=(1 + envelope)/2)
		lines(z, upper, lty=2, lwd=lwd, col=col.lines)
		lines(z, lower, lty=2, lwd=lwd, col=col.lines)
		if (line == "quartiles"){
			Q.x <- quantile(rstudent, c(.25,.75))
			Q.z <- if (distribution == 't') qt(c(.25,.75), df=res.df - 1) else qnorm(c(.25,.75))
			b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
			a <- Q.x[1] - b*Q.z[1]
			abline(a, b, col=col.lines, lwd=lwd)
		}
		if (line=="robust"){
			coef <- coefficients(rlm(ord.x~z))
			a <- coef[1]
			b <- coef[2]
			abline(a, b, col=col.lines, lwd=lwd)
		}                   
		result <- showLabels(z, ord.x,labels=ord.lab,  
				id.method = id.method, id.n = id.n, id.cex=id.cex, id.col=id.col)
	}
	if (length(result) == 0) invisible(result) else if (is.numeric(result)) sort(result) else result
}

qqPlot.glm <- function(x, ...){
	stop("QQ plot for studentized residuals not available for glm")
}
