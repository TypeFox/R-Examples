predcomp.test <- function(x, xhat1, xhat2, alternative = c("two.sided", "less", "greater"),
    lossfun = "losserr", lossfun.args = NULL, test = c("DM", "HG"), ...) {

    test <- match.arg( test )

    alternative <- tolower( alternative )
    alternative <- match.arg( alternative )

   out <- list()
   out$call <- match.call()

   e1 <- do.call(lossfun, c(list(x=x, xhat=xhat1), lossfun.args))
   e2 <- do.call(lossfun, c(list(x=x, xhat=xhat2), lossfun.args))

   if( test == "DM" ) {

	d <- e1 - e2
        d.cov.obj <- acf(d, type="covariance", plot=FALSE, na.action = na.pass, ...)
        d.cov <- d.cov.obj$acf[,,1]

	out$method <- "Diebold-Mariano Test"
	out$fitmodel <- "none"

	n <- length(d)
	d.var <- sum(c(d.cov[ 1 ], 2 * d.cov[-1])) / n

	STATISTIC <- mean(d, na.rm=TRUE)/sqrt(d.var)

	if (alternative == "two.sided") PVAL <- 2 * pnorm(-abs(STATISTIC))
   	else if (alternative == "less") PVAL <- pnorm(STATISTIC)
   	else if (alternative == "greater") PVAL <- pnorm(STATISTIC, lower.tail = FALSE)

   } else {

	fit <- hg.test( e1, e2, type = "OLS" )

	STATISTIC <- fit[ 1 ]
	PVAL <- fit[ 2 ]

	out$method <- "Hering-Genton Test"
	out$fitmodel <- "exponential"

   }

   alternative <- match.arg(alternative)

   out$loss.function <- lossfun
   out$loss.function.args <- lossfun.args
   out$statistic <- STATISTIC
   out$alternative <- alternative
   out$p.value <- PVAL
   out$data.name <- c(deparse(substitute(x)), deparse(substitute(xhat1)), deparse(substitute(xhat2)))

   class(out) <- c("predcomp.test", "htest")
   return(out)

} # end of 'predcomp.test' function.

losserr <- function(x, xhat, method=c("abserr","sqerr","simple","power","corrskill","dtw"),
    scale=1, p=1, dtw.interr=c("abserr","sqerr","simple","power"), ...) {

   method <- match.arg(method)
   if(method=="abserr") return(abs((xhat - x)/scale))
   else if(method=="sqerr") return(((xhat - x)/scale)^2)
   else if(method=="simple") return((xhat - x)/scale)
   else if(method=="power") return(((xhat - x)/scale)^p)
   else if(method=="corrskill") return(scale*(x - mean(x,na.rm=TRUE))*(xhat - mean(xhat, na.rm=TRUE)))
   else if(method=="dtw") {

	dtw.interr <- match.arg(dtw.interr)
	a <- dtw(xhat, x, step.pattern=asymmetric, ...)
	w <- numeric(max(a$index1, a$index2)) + NA
	w[a$index2] <- xhat[a$index1]
	d1 <- abs(a$index1 - a$index2)
	if(dtw.interr=="abserr") d2 <- abs((w - x)/scale)
	else if(dtw.interr=="sqerr") d2 <- ((w - x)/scale)^2
	else if(dtw.interr=="simple") d2 <- (w - x)/scale
	else if(dtw.interr=="power") d2 <- ((w - x)/scale)^p
	else stop("losserr: dtw.interr must be one of abserr, sqerr, simple, or power")

	return(d1 + d2)

   } else stop("losserr: method must be one of abserr, sqerr, simple, power, corrskill, or dtw")

} # end of 'losserr' function.

exponentialACV <- function(x, y, ...) {

   args <- list(...)
   if(!is.null(args$start)) res <- nls(y~sigma^2*exp(-3*x/theta), data=data.frame(x=x, y=y), ...)
   else {

	if(any(y<0.05)) theta.start <- min(x[y<0.5],na.rm=TRUE)
	else theta.start <- length(x)
	res <- nls(y~sigma^2*exp(-3*x/theta), data=data.frame(x=x, y=y), start=list(sigma=sqrt(y[1]), theta=theta.start), ...)

   }

   return(res)

} # end of 'acvparametric' function.

summary.predcomp.test <- function(object, ...) {

   a <- object

   print(a$call)

   cat("\n", "Loss function used is: ", a$loss.function, "\n")

   if(!is.null(a$loss.function.args)) if(!is.null(a$loss.function.args$method)) {

	m <- a$loss.function.args$method
	if(m == "abserr") msg <- "Absolute Error Loss"
	else if(m == "sqerr") msg <- "Square Error Loss"
	else if(m == "simple") msg <- "Simple Error Loss"
	else if(m == "power") {
	   if(is.null(a$loss.function.args$p)) p <- 1
	   else p <- a$loss.function.args$p
	   msg <- paste("Power Error loss with p = ", p, sep="")

	} else if(m == "corrskill") msg <- "Correlation Skill"
	else if(m == "dtw") {

	   if(is.null(a$loss.function.args$dtw.interr)) interr <- "Absolute Error Loss"
	   else if(a$loss.function.args$dtw.interr == "abserr") interr <- "Absolute Error Loss"
	   else if(a$loss.function.args$dtw.interr == "sqerr") interr <- "Square Error Loss"
	   else if(a$loss.function.args$dtw.interr == "simple") interr <- "Simple Error Loss"
	   else if(a$loss.function.args$dtw.interr == "power") {
		if(is.null(a$loss.function.args$p)) p <- 1
           	else p <- a$loss.function.args$p
           	msg <- paste("Power Error loss with p = ", p, sep="")

	   }

	   msg <- paste("Discrete Time Warping Loss with Intensity Error = ", interr, sep="")

	}

	cat(paste("Method used: ", msg, sep=""), "\n")
	a$loss.message <- msg

   }

   if(is.null(a$loss.function.args)) if(a$loss.function == "losserr") cat("Absolute Error Loss\n")

   print(a)

   invisible(a)

} # end of 'summary.predcomp.test' function.
