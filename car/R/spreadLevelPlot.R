# spread-level plots (J. Fox)

# 16 March 2010 by J. Fox: spreadLevelPlot.lm now deletes observations with negative fitted values
# 25 May 2010 by J. Fox: corrected errors due to introduction of grid()
# 2015-11-24: added smoother and related args to lm method. John

slp <- function(...) spreadLevelPlot(...)

spreadLevelPlot <- function(x, ...) {
	UseMethod("spreadLevelPlot")
}

spreadLevelPlot.default <- function(x, by, robust.line=TRUE, 
	start=0, xlab="Median", ylab="Hinge-Spread", point.labels=TRUE, las=par("las"),
	main=paste("Spread-Level Plot for", deparse(substitute(x)), 
		"by", deparse(substitute(by))), col=palette()[1], col.lines=palette()[2],
    pch=1, lwd=2, grid=TRUE, ...){
	good <- complete.cases(x, by)
	if (sum(good) != length(x)) {
		warning("NAs ignored")
		x <- x[good]
		by <- by[good]
	}    
	min.x <- min(x)
	if (min.x <= -start){
		start <- nice(-min.x + 0.05*diff(quantile(x, c(.25, .75))), direction="up")
		warning(paste("Start =", start," added to avoid 0 or negative values."))
	}
	if (start != 0) {
		xlab <- paste(xlab, "+", signif(start, getOption("digits")))
		x <- x + start
	}
	values <- unique(as.character(by))
	result <- matrix(0, length(values), 4)
	dimnames(result) <-list(values, c("LowerHinge", "Median", "UpperHinge", "Hinge-Spread"))
	for (i in seq(along=values)){
		five <- fivenum(x[by == values[i]])
		result[i, ] <- c(five[2:4], five[4] - five[2])
	}
	medians<-result[ ,2]
	spreads<-result[ ,4]
	plot(medians, spreads, type="n", log="xy", main=main, xlab=xlab, ylab=ylab, 
		las=las, pch=pch, col=col, ...)
	if(grid){
    grid(lty=1, equilogs=FALSE)
    box()}
	points(medians, spreads, col=col, pch=pch)
	pos <- ifelse(medians > median(medians), 2, 4)
	if (point.labels) text(medians, spreads, as.character(values), pos=pos, ...)
	mod <- if (robust.line)
		rlm(log(spreads) ~ log(medians))
	else lm(log(spreads) ~ log(medians), ...)
	ord <- order(medians)
	first <- ord[1]
	last <- ord[length(ord)]
	lines(start + medians[c(first, last)], exp(fitted.values(mod)[c(first, last)]), 
		col=col.lines, lwd=lwd, ...)
	p <- 1 - (coefficients(mod))[2]
	names(p) <- NULL
	result <- list(Statistics=as.data.frame(result[ord,]), PowerTransformation=p)
	class(result) <- "spreadLevelPlot"
	result
}

spreadLevelPlot.lm <- function(x, robust.line=TRUE, 
    smoother=loessLine, smoother.args=list(),
	xlab="Fitted Values",
	ylab="Absolute Studentized Residuals", las=par("las"),
	main=paste("Spread-Level Plot for\n", deparse(substitute(x))),
	pch=1, col=palette()[1], col.lines=palette()[2], col.smoother=palette()[3], lwd=2, grid=TRUE, ...){
	resid <- na.omit(abs(rstudent(x)))
	fitval <- na.omit(fitted.values(x))
	non.pos <- fitval <= 0
	if (any(non.pos)){
		fitval <- fitval[!non.pos]
		resid <- resid[!non.pos]
		n.non.pos <- sum(non.pos)
		warning(n.non.pos, " negative", if(n.non.pos > 1) " fitted values" else " fitted value", " removed")
	}
	min <- min(fitval)
	plot(fitval, resid, log="xy", main=main, xlab=xlab, ylab=ylab, 
			las=las, col=col, pch=pch, type="n", ...)
	if(grid){
    grid(lty=1, equilogs=FALSE)
    box()}
	points(fitval, resid, col=col, pch=pch)
	mod <- if (robust.line)
			rlm(log(resid) ~ log(fitval))
		else lm(log(resid) ~ log(fitval), ...)
	first <- which.min(fitval) 
	last <- which.max(fitval) 
	lines((fitval)[c(first, last)], exp(fitted.values(mod)[c(first, last)]), 
		lwd=lwd, col=col.lines, ...)
	if (is.null(smoother.args$lwd)) smoother.args$lwd <- lwd
	if (is.null(smoother.args$lty)) smoother.args$lty <- 2
	if (is.function(smoother)) smoother(fitval, resid, col=col.smoother,
	    log.x=TRUE, log.y=TRUE, smoother.args=smoother.args)
	p <- 1 - (coefficients(mod))[2]
	names(p) <- NULL
	result <- list(PowerTransformation=p)
	class(result) <- "spreadLevelPlot"
	result
}

spreadLevelPlot.formula <- function (x, data=NULL, subset, na.action, 
	main=paste("Spread-Level Plot for", varnames[response], "by", varnames[-response]), ...) {
	if (missing(na.action)) 
		na.action <- getOption("na.action")
	m <- match.call(expand.dots = FALSE)
	m$formula <- x
	if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
		m$data <- as.data.frame(data)
	m$... <- m$main <- m$x <- NULL
	m[[1]] <- as.name("model.frame")
	mf <- eval(m, sys.frame(sys.parent()))
	response <- attr(attr(mf, "terms"), "response")
	varnames <- names(mf)
	if (!response) stop ("no response variable specified")
	if (length(varnames) > 2) stop("right-hand side of model has more than one variable")
	x <- mf[[response]]
	by <- mf[[varnames[-response]]]
	spreadLevelPlot(x, by, main=main, ...)
}

print.spreadLevelPlot <- function(x, ...){
	if (!is.null(x$Statistics)) print(x$Statistics, ...)
	cat('\nSuggested power transformation: ', x$PowerTransformation,'\n')
	invisible(x)
}
