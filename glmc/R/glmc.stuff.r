print.summary.glmc <-
    function (x, digits = max(3, getOption("digits") - 3),
	      symbolic.cor = x$symbolic.cor,
	      signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")
    cat("Deviance Residuals: \n")
    if(x$df.residual > 5) {
	x$deviance.resid <- quantile(x$deviance.resid,na.rm=TRUE)
	names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
    }
    print.default(x$deviance.resid, digits=digits, na.print = "", print.gap = 2)

    if(length(x$aliased) == 0) {
        cat("\nNo Coefficients\n")
    } else {
        ## df component added in 1.8.0
        if (!is.null(df<- x$df) && (nsingular <- df[3] - df[1]))
            cat("\nCoefficients: (", nsingular,
                " not defined because of singularities)\n", sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if(!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4,
                            dimnames=list(cn, colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                     na.print="NA", ...)
    }
    ##
    cat("\n(Dispersion parameter for ", x$family$family,
	" family taken to be ", format(x$dispersion), ")\n\n",
	apply(cbind(paste(format.default(c("Null","Residual"),width=8,flag=""),
			  "deviance:"),
		    format(unlist(x[c("null.deviance","deviance")]),
			   digits= max(5, digits+1)), " on",
		    format(unlist(x[c("df.null","df.residual")])),
		    " degrees of freedom\n"),
	      1, paste, collapse=" "),
	"AIC: ", format(x$aic, digits= max(4, digits+1)),"\n\n",
	"Number of Fisher Scoring iterations: ", x$iter,
	"\n", sep="")

    correl <- x$correlation
    if(!is.null(correl)) {
# looks most sensible not to give NAs for undefined coefficients
#         if(!is.null(aliased) && any(aliased)) {
#             nc <- length(aliased)
#             correl <- matrix(NA, nc, nc, dimnames = list(cn, cn))
#             correl[!aliased, !aliased] <- x$correl
#         }
	p <- NCOL(correl)
	if(p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
		print(symnum(correl, abbr.colnames = NULL))
	    } else {
		correl <- format(round(correl, 2), nsmall = 2, digits = digits)
		correl[!lower.tri(correl)] <- ""
		print(correl[-1, -p, drop=FALSE], quote = FALSE)
	    }
	}
    }
    cat("\n")
    invisible(x)
}


## GLM Methods for Generic Functions :

coef.glmc <- function(object, ...) object$coefficients
deviance.glmc <- function(object, ...) object$deviance
effects.glmc <- function(object, ...) object$effects
fitted.glmc <- function(object, ...)
{
    if(is.null(object$na.action)) object$fitted.values
    else napredict(object$na.action, object$fitted.values)
}

family.glmc <- function(object, ...) object$family

residuals.glmc <-
    function(object,
	     type = c("deviance", "pearson", "working", "response", "partial"),
	     ...)
{
    type <- match.arg(type)
    y <- object$y
    r <- object$residuals
    mu	<- object$fitted.values
    wts <- object$prior.weights
    res <- switch(type,
		  deviance = if(object$df.res > 0) {
		      d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, wts), 0))
		      ifelse(y > mu, d.res, -d.res)
		  } else rep.int(0, length(mu)),
		  pearson = (y-mu)*sqrt(wts)/sqrt(object$family$variance(mu)),
		  working = r,
		  response = y - mu,
		  partial = r + predict(object,type="terms")
		  )
    if(is.null(object$na.action)) res
    else naresid(object$na.action, res)
}
