quantilefun <-
    function(y)
    approxfun(seq(0, 1, length = length(y)), sort(y),
              yleft = NA, yright = NA)

percentilefun <-
    function(y)
    approxfun(sort(y), seq(0, 1, length = length(y)),
              yleft = 0, yright = 1)

qqnorm <-
    function(y, pch = 20,
             xlab = "Standard Normal Quantiles",
             ylab = "Sample Quantiles", make.plot=TRUE, ...)
    {
	args <- list(...)
        y <- sort(na.omit(y))
        n <- length(y)
        p <- (1:length(y) - .5)/length(y)
        k <- .895 / (sqrt(n) * (1 - .01 / sqrt(n) + .85 / n))
        l <- suppressWarnings(qnorm(p - k))
        q <- qnorm(p)
        u <- suppressWarnings(qnorm(p + k))
	if(make.plot) {
           if(is.null(args$xlim)) plot(q, y, xlim = range(l, q, u, na.rm = TRUE), xlab = xlab, ylab = ylab, pch = pch, ...)
	   else plot(q, y, xlab = xlab, ylab = ylab, pch = pch, ...)
           lines(l, y, lty = 2, col = "darkgray")
           lines(u, y, lty = 2, col = "darkgray")
	}
	out <- data.frame(lower=l, upper=u, qnorm=q, data=y)
	invisible(out)
    }

qqplot <-
    function(x, y, pch = 20,
             xlab = "x Quantiles",
             ylab = "y Quantiles", regress = TRUE, make.plot=TRUE, ...)
    {
	args <- list(...)
	out <- list()
	out$call <- match.call()
	out$names <- list(x = as.character(substitute(x)), y = as.character(substitute(y)))
        x <- sort(na.omit(x))
        y <- sort(na.omit(y))
        qy <- quantilefun(y)
        m <- length(x)
        n <- length(y)
        N <- m + n
        M <- m * (n / N)
        K <- 1.36
        p <- (1:m - 1)/(m - 1)
        yq <- qy(p)
        yl <- qy(p - K/sqrt(M))
        yu <- qy(p + K/sqrt(M))
	if(make.plot) {
           if(is.null(args$xlim) && is.null(args$ylim)) plot(x, yq, pch = pch, xlim = range(x), ylim = range(yq, yl, yu, na.rm = TRUE), xlab = xlab, ylab = ylab, ...)
	   else if(is.null(args$xlim)) plot(x, yq, pch = pch, xlim = range(x), xlab = xlab, ylab = ylab, ...)
	   else if(is.null(args$ylim)) plot(x, yq, pch = pch, ylim = range(yq, yl, yu, na.rm = TRUE), xlab = xlab, ylab = ylab, ...)
	   else plot(x, yq, pch = pch, xlab = xlab, ylab = ylab, ...)
           lines(x, yl, lty = 2, col = "gray")
           lines(x, yu, lty = 2, col = "gray")
	   abline(0,1, lty=2, col="darkorange")
	}
	if(regress) {
	   fit <- lm(y~x, data=data.frame(x=x, y=yq))
	   if(make.plot) {
		lines(x, predict(fit), col="grey", lty=1)
		legend("topleft", legend=c("1-1 line", "regression line", "95% confidence bands"), col=c("darkorange","grey","gray"), lty=c(2,1,2), bty="n")
	   }
	   out$regression <- fit
	} else if(make.plot) legend("bottomright", legend=c("1-1 line", "95% confidence bands"), col=c("darkorange","gray"), lty=c(2,2), bty="n")
	out$qdata <- data.frame(x=x, y=yq, lower=yl, upper=yu)
	class(out) <- "qqplot"
	invisible(out)
    } # end of 'qqplot' function.

summary.qqplot <- function(object, ...) {
   print(object$call)
   if(!is.null(object$regression)) print(summary(object$regression))
   invisible(object)
} # end of 'summary.qqplot' function.

plot.qqplot <- function(x, ...) {
   args <- list(...)
   if(is.null(args$xlab)) xlab <- "x Quantiles"
   if(is.null(args$ylab)) ylab <- "y Quantiles"
   z <- x$qdata$x
   yq <- x$qdata$y
   yl <- x$qdata$lower
   yu <- x$qdata$upper
   if(is.null(args$xlab) && is.null(args$ylab)) {
      if(is.null(args$xlim) && is.null(args$ylim)) plot(z, yq, xlim = range(z), ylim = range(yq, yl, yu, na.rm = TRUE), xlab = xlab, ylab = ylab, ...)
      else if(is.null(args$xlim)) plot(z, yq, xlim = range(z), xlab = xlab, ylab = ylab, ...)
      else if(is.null(args$ylim)) plot(z, yq, ylim = range(yq, yl, yu, na.rm = TRUE), xlab = xlab, ylab = ylab, ...)
      else plot(z, yq, xlab = xlab, ylab = ylab, ...)
   } else if(is.null(args$xlab)) {
      if(is.null(args$xlim) && is.null(args$ylim)) plot(z, yq, xlim = range(z), ylim = range(yq, yl, yu, na.rm = TRUE), xlab = xlab, ...)
      else if(is.null(args$xlim)) plot(z, yq, xlim = range(z), xlab = xlab, ...)
      else if(is.null(args$ylim)) plot(z, yq, ylim = range(yq, yl, yu, na.rm = TRUE), xlab = xlab, ...)
      else plot(z, yq, xlab = xlab, ...)
   } else if(is.null(args$ylab)) {
      if(is.null(args$xlim) && is.null(args$ylim)) plot(z, yq, xlim = range(z), ylim = range(yq, yl, yu, na.rm = TRUE), ylab = ylab, ...)
      else if(is.null(args$xlim)) plot(z, yq, xlim = range(z), ylab = ylab, ...)
      else if(is.null(args$ylim)) plot(z, yq, ylim = range(yq, yl, yu, na.rm = TRUE), ylab = ylab, ...)
      else plot(z, yq, ylab = ylab, ...)
   } else {
      if(is.null(args$xlim) && is.null(args$ylim)) plot(z, yq, xlim = range(z), ylim = range(yq, yl, yu, na.rm = TRUE), ...)
      else if(is.null(args$xlim)) plot(z, yq, xlim = range(z), ...)
      else if(is.null(args$ylim)) plot(z, yq, ylim = range(yq, yl, yu, na.rm = TRUE), ...)
      else plot(z, yq, ...)
   }
   lines(z, yl, lty = 2, col = "gray")
   lines(z, yu, lty = 2, col = "gray")
   abline(0,1, lty=2, col="darkorange")
   if(!is.null(x$regression)) {
	fit <- x$regression
	lines(z, predict(fit), col="grey", lty=1, lwd=1.5)
        legend("bottomright", legend=c("1-1 line", "regression line", "95% confidence bands"), col=c("darkorange","grey","gray"), lty=c(2,1,2), bty="n")
   } else legend("bottomright", legend=c("1-1 line", "95% confidence bands"), col="darkorange", lty=c(2,1), bty="n")
   invisible(x)
} # end of 'plot.qqplot' function.

shiftplot <- function(x, y, pch = 20, xlab = "x Quantiles", ylab = "y Quantiles", main = NULL, ...) {
        x <- sort(na.omit(x))
        y <- sort(na.omit(y))
        qy <- quantilefun(y)
        m <- length(x)
        n <- length(y)
        N <- m + n
        M <- m * n / N
        K <- 1.36
        p <- (1:m - 1)/(m - 1)
        yq <- qy(p) - x
        yl <- qy(p - K/sqrt(M)) - x
        yu <- qy(p + K/sqrt(M)) - x
        plot(x, yq, pch = pch, xlim = range(x), ylim = range(yq, yl, yu, na.rm = TRUE), xlab = xlab, ylab = ylab, main = main, ...)
        lines(x, yl, col = "darkgray")
        lines(x, yu, col = "darkgray")
    } # end of 'shiftplot' function.

# qqevd <- function( object, conf=NULL, ...) {
#    if( class( object) == "gev.fit") {
# 	n <- length( object$data)
# 	x <- (1:n)/(n + 1)
# 	if( object$trans) {
# 	   plot( -log(-log(x)), sort( object$data), xlab="Empirical", ylab="Model",
# 			main="Residual Quantile Plot (Gumbel Scale)", ...)
# 	   abline(0,1,col=4,lwd=1.5)
# 	   if( !is.null( conf)) {
# 		## Need to make a change in ismev to allow one to know which are the parameter values
# 		## for the case where all covariates are zero!  Then, change mle below accordingly.
# 		d <- gev.rl.gradient( a=object$mle, p=1-1/(log(x)))
# 		v <- apply( d, 1, q.form, m=object$cov)
# 	   }
# 	} # end of if 'trans' stmts.
#    } else if( class( object) == "gpd.fit") {
# 
#    } else {
# 
#    }
# } # end of 'qqevd' function.
