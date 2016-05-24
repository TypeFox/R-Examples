##'Standard accessor functions for ADMB model fits
##'
##'Extract standard information such as log-likelihood, AIC, coefficients, etc.
##'from ADMB model fits
##'

##' @export

##'@param x an ADMB model fit (of class "admb")
##'@param object an ADMB model fit (of class "admb")
##'@param k penalty value for AIC fits
##'@param type which type of parameters to report. Character vector, including
##' one or more of "fixed" or "par" (standard, fixed-effect parameters);
##' "random" (random effect parameters); "rep" (report variables); "sdrpt" (sdreport variables);
##' "extra" (report and sdreport); "all" (all of the above).
##'@param parm (currently ignored: FIXME) select parameters
##'@param level alpha level for confidence interval
##'@param method (character): "default" or "quad", quadratic (Wald) intervals
##'based on approximate standard errors; "profile", profile CIs (if profile was
##'computed); "quantile", CIs based on quantiles of the MCMC-generated posterior
##'density (if MCMC was computed); "HPDinterval", CIs based on highest posterior
##'density (ditto)
##'@param correlation currently unused parameter
##'@param symbolic.cor currently unused parameter
##'@param verbose show messages
##'@param \dots other parameters (for S3 generic compatibility)
##'@return Extracts appropriate values: numeric (scalar) for AIC, type logLik
##'for logLik, numeric vector of coefficients, numeric variance-covariance
##'matrix of parameter estimates
##'@author Ben Bolker
##'@keywords misc
##'@examples
##'
##'  admbex <- system.file("doc","Reedfrog_runs.RData",package="R2admb")
##'  load(admbex)
##'  m1
##'  coef(m1)
##'  summary(m1)
##'  coef(summary(m1)) ## returns just z-table
##'  AIC(m1)
##'  vcov(m1)
##'  logLik(m1)
##'  deviance(m1)
##'  stdEr(m1)
##'
AIC.admb <- function(object,...,k=2) {
	if (length(list(...))>0) stop("multi-object AIC not yet implemented")
	deviance(object)+k*length(coef(object))
}

## copied from stats::
format.perc <- function (probs, digits)  {
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), 
          "%")
}

##' @rdname AIC.admb
##' @export
##' @importFrom coda HPDinterval as.mcmc

confint.admb <- function(object, parm, level=0.95, method="default", type="fixed", ...) {
    if (!missing(parm) && is.character(parm)) {
        ## try to catch mistakes like specifying method or type in parm slot
        if (any(is.na(m <- match(parm,names(coef(object,"all")))))) {
            missparms <- parm[is.na(m)]
            stop("requested parameters missing from parameter vector: ",
                 paste(missparms,collapse=", "))
        }
    }
    if (method %in% c("default","quad")) {
        ## copied from confint.default because we want to keep the *default* type
        ## for vcov() equal to "fixed", and can't pass options through confint.default()
        ## as long as we're at it we should use $se rather than sqrt(diag(vcov))
        cf <- coef(object,type="all")
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        pct <- format.perc(a, 3)
        fac <- qnorm(a)
        ci <- array(NA, dim = c(length(cf), 2L), dimnames = list(names(cf),pct))
        ses <- object$se
        ci[] <- cf + ses %o% fac
        tab <- ci[get_parn(object,type),]
    } else if (method=="profile") {
        vals <- object[["prof"]]
        if (is.null(vals)) stop("model not fitted with profile=TRUE")
        if (!level %in% c(0.9,0.95,975)) stop("arbitrary levels not yet implemented:",
                                              "level must be in (0.9,0.95,0.975)")
        tab <- t(sapply(vals,function(x) {
            x$ci[x$ci[,"sig"]==level,c("lower","upper")]
        }))
        colnames(tab) <- paste(c((1-level)/2,(1+level)/2)*100,"%")
    } else if (method %in% c("quantile","HPDinterval")) {
        vals <- object[["mcmc"]]
        if (is.null(vals)) stop("model not fitted with mcmc=TRUE")
        if (method=="quantile") {
            tab <- t(apply(vals,2,quantile,c((1-level)/2,(1+level)/2)))
        } else {
            tab <- HPDinterval(as.mcmc(vals))
            colnames(tab) <- paste(c((1-level)/2,(1+level)/2)*100,"%")
        }
    }
    if (missing(parm)) parm <- seq(nrow(tab))
    tab[parm,,drop=FALSE]
}

##' @rdname AIC.admb
##' @export
print.admb <- function(x, verbose=FALSE, ...) {
	cat("Model file:",x$fn,"\n")
	if (is.null(x$loglik)) {
		cat("No fit\n")
		return(invisible(NULL))
	}
	cat("Negative log-likelihood:",-x$loglik,"\n")
	cat("Coefficients:\n")
	print(coef(x))
	## FIXME: indicate extra parameters?
	if (!is.null(x$mcmc)) {
		mcpar <- attr(x$mcmc,"mcpar")
		cat("MCMC parameters: start=",mcpar[1],", end=",mcpar[2],", thin=",mcpar[3],"\n",sep="")
	}
	if (verbose) cat(x$txt,sep="\n")
}      

##' @rdname AIC.admb
##' @export

summary.admb <- function(object, correlation=FALSE, symbolic.cor = FALSE, ...) {
        coef.p <- unlist(coef(object,"par"))
	s.err <- sqrt(diag(vcov(object)))
	tvalue <- coef.p/s.err
	dn <- c("Estimate", "Std. Error")
	pvalue <- 2 * pnorm(-abs(tvalue))
	coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
	dimnames(coef.table) <- list(names(coef.p), c(dn, 
					"z value", "Pr(>|z|)"))
	ans <- c(list(coefficients=coef.table),
                 object[c("loglik","fn","npar")])
	class(ans) <- "summary.admb"
	ans
}

##' @rdname AIC.admb
##' @export
##' @param digits number of digits to display
##' @param signif.stars show significance stars?
print.summary.admb <- function(x,
		digits = max(3, getOption("digits") - 3),
		symbolic.cor = x$symbolic.cor, 
		signif.stars = getOption("show.signif.stars"), ...) {
	coefs <- x$coefficients
	cat("Model file:",x$fn,"\n")
	cat("Negative log-likelihood: ",sprintf("%1.1f",-x$loglik),"\t",
            "AIC: ",sprintf("%.1f",-2*(x$loglik-x$npar)),"\n")
	cat("Coefficients:\n")
	printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
			na.print = "NA", ...)
}

## utility function: get positions in parameter vector matching various types
get_parn <- function(x,type=c("par","fixed","random","extra","sdrpt","rep","all")) {
    type <- match.arg(type)
    ## "par" and "fixed": synonyms
    if (type=="all") return(seq(x$npar_total))
    type[type=="par"] <- "fixed"
    w <- numeric(0)
    if ("extra" %in% type && any(c("sdrpt","rep") %in% type))
        stop("both 'extra' and 'sdrpt/rep' specified")
    sseq <- function(n) if (n==0) numeric(0) else seq(n)
    if ("fixed" %in% type) w <-c(w,sseq(x$npar))
    if ("random" %in% type) w <- c(w,x$npar+sseq(x$npar_re))
    if ("sdrpt" %in% type) w <- c(w,x$npar+x$npar_re+sseq(x$npar_sdrpt))
    if ("rep" %in% type) w <- c(w,x$npar+x$npar_re+x$npar_sdrpt+sseq(x$npar_rep))
    if ("extra" %in% type) w <- c(w,x$npar+x$npar_re+sseq(x$npar_sdrpt+x$npar_rep))
    w
}
    
##' @rdname AIC.admb
##' @export
logLik.admb <- function(object,...) {
    L <- object$loglik
    df <- length(coef(object))
    attr(L,"df") <- df
    class(L) <- "logLik"
    ## fixme: would be nice to have an "nobs" attribute
    ##   but not sure when/how it can be defined ...
    L
}

##' @rdname AIC.admb
##' @export
coef.admb <- function(object,type="fixed",...) {
    object$coefficients[get_parn(object,type)]
}

##' @rdname AIC.admb
##' @export
vcov.admb <- function(object,type="fixed",...) {
    v <- object$vcov
    w <- get_parn(object,type)
    v[w,w,drop=FALSE]
}

##' @rdname AIC.admb
##' @export stdEr
stdEr <- function(object, ...) {
    UseMethod("stdEr")
}

##' @rdname AIC.admb
##' @export
stdEr.admb <- function(object,type="fixed",...) {
    s <- sqrt(diag(object$vcov))
    object$se[get_parn(object,type)]
}

##' @rdname AIC.admb
##' @export
deviance.admb <- function(object,...) -2*object$loglik
