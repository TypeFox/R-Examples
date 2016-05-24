################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Methods for "twinSIR" fits, specifically:
### - vcov: enabling the use of function confint to calculate Wald
###         confidence intervals for the parameter estimates.
### - logLik: enables the use of function AIC
### - AIC, extractAIC: compute AIC or OSAIC depending on argument 'one.sided'
### - print, summary, print.summary, plot (intensityPlot), ...
###
### Copyright (C) 2009-2014 Sebastian Meyer, contributions by Michael Hoehle
### $Revision: 1088 $
### $Date: 2014-10-24 09:29:43 +0200 (Fre, 24. Okt 2014) $
################################################################################

### don't need a specific coef-method (identical to stats:::coef.default)
## coef.twinSIR <- function (object, ...)
## {
##     object$coefficients
## }

# asymptotic variance-covariance matrix (inverse of fisher information matrix)
vcov.twinSIR <- function (object, ...)
{
    solve(object$fisherinfo)
}

logLik.twinSIR <- function (object, ...)
{
    r <- object$loglik
    attr(r, "df") <- length(coef(object))
    class(r) <- "logLik"
    r
}

# Note: pz is determined by scanning the names of coef(object),
#       thus the 'model' component is not necessary
# See the Hughes and King (2003) paper for details
.OSAICpenalty <- function (twinSIRobject, k = 2, nsim = 1e3)
{
    theta <- coef(twinSIRobject)
    npar <- length(theta)
    pz <- length(grep("cox\\([^)]+\\)", names(theta), ignore.case = FALSE,
                      perl = FALSE, fixed = FALSE, useBytes = FALSE,
                      invert = FALSE))
    px <- npar - pz   # number of constrained (non-negative) parameters
    
    penalty <- if (px == 0L) {
        k * pz   # default AIC penalty (with k = 2)
    } else if (px == 1L) {
        k * (pz + 0.5)
    } else if (px == 2L) {
        Sigma <- vcov(twinSIRobject)   # parameter covariance matrix
        rho <- cov2cor(Sigma[1:2,1:2])[1,2]
        as <- acos(rho)/2/pi
        w <- c(as, 0.5, 0.5-as)
        k * sum(w * (pz + 0:2))   # = k * sum(w * (npar - px + 0:2))
    } else { # px > 2
        message("Computing OSAIC weights for ", px,
                " epidemic covariates based on ", nsim, " simulations ...")
        W <- vcov(twinSIRobject)[1:px,1:px]
        w.sim <- w.chibarsq.sim(p=px, W=W, N=nsim)
          #c.f. (12) in Hughes & King (2003), r_i=px, m=0:px, ki=npar
          #as npar=pz+px, we have that npar-px = pz, hence the sum is
        k * sum(w.sim * (pz + 0:px))
    }
    
    attr(penalty, "exact") <- px <= 2
    penalty
}

AIC.twinSIR <- function (object, ..., k = 2, one.sided = NULL, nsim = 1e3)
{
    AIC.default <- match.call()
    AIC.default$one.sided <- NULL
    AIC.default$nsim <- NULL
    AIC.default[[1]] <- call(":::", as.name("stats"), as.name("AIC.default"))
    ## I don't see any easy way of using AIC.default while avoiding ":::".
    ## NextMethod() does not fit due to extra arguments one.sided & nsim.
    ## Could maybe unclass "object" and all objects in "..." and then use AIC()
    
    if (is.null(one.sided)) {
        one.sided <- object$method == "L-BFGS-B"
    }
    
    if (one.sided) {
        penalty <- .OSAICpenalty(object, k = k, nsim = nsim)
        edf <- length(coef(object))
        AIC.default$k <- penalty/edf
    }

    res <- eval(AIC.default, parent.frame())
    attr(res, "type") <- if (one.sided) "One-sided AIC" else "Standard AIC"
    attr(res, "exact") <- if (one.sided) attr(penalty, "exact") else TRUE
    res
}

extractAIC.twinSIR <- function (fit, scale = 0, k = 2, one.sided = NULL,
    nsim = 1e3, ...)
{
    if (is.null(one.sided)) {
        one.sided <- fit$method == "L-BFGS-B"
    }
    
    loglik <- logLik(fit)
    edf <- attr(loglik, "df")
    penalty <- if (one.sided) {
                   .OSAICpenalty(fit, k = k, nsim = nsim)   # one-sided AIC
               } else {
                   k * edf                                  # default AIC
               }
    res <- c(edf = edf, AIC = -2 * c(loglik) + penalty)            
    
    attr(res, "type") <- if (one.sided) "One-sided AIC" else "Standard AIC"
    attr(res, "exact") <- if (one.sided) attr(penalty, "exact") else TRUE
    res
}

print.twinSIR <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")
    print.default(x$call)
    cat("\nCoefficients:\n")
    print.default(format(coef(x), digits=digits), print.gap = 2, quote = FALSE)
    cat("\nLog-likelihood: ", format(logLik(x), digits=digits), "\n", sep = "")
    if (!x$converged) {
        cat("\nWARNING: OPTIMIZATION DID NOT CONVERGE!\n")
    }
    cat("\n")
    invisible(x)
}

summary.twinSIR <- function (object,
    correlation = FALSE, symbolic.cor = FALSE, ...)
{
    ans <- object[c("call", "converged", "counts", "intervals", "nEvents")]
    ans$cov <- vcov(object)
    est <- coef(object)
    se <- sqrt(diag(ans$cov))
    zval <- est/se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    ans$coefficients <- cbind(est, se, zval, pval)
    dimnames(ans$coefficients) <- list(names(est),
        c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    if (correlation) {
        ans$correlation <- cov2cor(ans$cov)
        ans$symbolic.cor <- symbolic.cor
    }
    ans$loglik <- logLik(object)
    aic <- extractAIC(object, ...)
    ans$aic <- as.vector(aic[2L])   # remove 'edf' element
    attributes(ans$aic) <- attributes(aic)[c("type", "exact")]
    class(ans) <- "summary.twinSIR"
    ans
}

print.summary.twinSIR <- function (x,
    digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
    signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n")
    print.default(x$call)
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
    nEvents <- x$nEvents
    nh0 <- length(nEvents)
    if (nh0 < 2L) {
        cat("\nTotal number of infections: ", nEvents, "\n")
    } else {
        cat("\nBaseline intervals:\n")
        intervals <- character(nh0)
        for(i in seq_len(nh0)) {
            intervals[i] <-
            paste("(",
                  paste(format(x$intervals[c(i,i+1L)],trim=TRUE), collapse=";"),
                  "]", sep = "")
        }
        names(intervals) <- paste("logbaseline", seq_len(nh0), sep=".")
        print.default(rbind("Time interval" = intervals,
                            "Number of events" = nEvents),
                      quote = FALSE, print.gap = 2)
    }
    cat("\n", attr(x$aic, "type"), ": ", format(x$aic, digits=max(4, digits+1)),
        if (!attr(x$aic, "exact")) "\t(simulated penalty weights)" else "",
        sep = "")
    cat("\nLog-likelihood:", format(x$loglik, digits = digits))
    cat("\nNumber of log-likelihood evaluations:", x$counts[1], "\n")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1L) {
        cat("\nCorrelation of Coefficients:\n")
        if (is.logical(symbolic.cor) && symbolic.cor) {
            correl <- symnum(correl, abbr.colnames = NULL)
            correlcodes <- attr(correl, "legend")
            attr(correl, "legend") <- NULL
            print(correl)
            cat("---\nCorr. codes:  ", correlcodes, "\n", sep="")
        } else {
            correl <- format(round(correl, 2), nsmall = 2, digits = digits)
            correl[!lower.tri(correl)] <- ""
            print(correl[-1, -p, drop = FALSE], quote = FALSE)
        }
        }
    }
    if (!x$converged) {
        cat("\nWARNING: OPTIMIZATION DID NOT CONVERGE!\n")
    }
    cat("\n")
    invisible(x)
}


### Plot method for twinSIR (wrapper for intensityplot)

plot.twinSIR <- function (x, which, ...) # defaults for 'which' are set below
{
    cl <- match.call()
    cl[[1]] <- as.name("intensityplot")
    eval(cl, envir = parent.frame())
}

formals(plot.twinSIR)$which <- formals(intensityplot.twinSIR)$which



######################################################################
# Extract the "residual process" (cf. Ogata, 1988), i.e. the
# fitted cumulative intensity at the event times.
# -> "generalized residuals similar to those discussed in Cox and Snell (1968)"
######################################################################

residuals.twinSIR <- function(object, ...)
{
  #Extract event and stop-times
  eventTimes <- attr(object$model$survs,"eventTimes")
  sortedStop <- sort(unique(object$model$survs[,"stop"]))
  eventTimesIdx <- match(eventTimes, sortedStop)
  
  #Dimensions and zero vector (in case we need it)
  nTimes <- nrow(object$model$X)
  zerovec <- numeric(nTimes)

  # Extract the fitted model params
  px <- ncol(object$model$X)
  pz <- ncol(object$model$Z)
  theta <- coef(object)
  alpha <- theta[seq_len(px)]
  beta <- theta[px+seq_len(pz)]

  # Initialize e, h and thus lambda
  if (px > 0) { e <- as.vector(object$model$X %*% as.matrix(alpha)) } else { e <- zerovec }
  if (pz > 0) { h <- as.vector(exp(object$model$Z %*% as.matrix(beta))) } else { h <- zerovec }
  lambda <- (e + h)

  #Determine bloks
  BLOCK <- as.numeric(factor(object$model$survs$start))

  # lambda_i integrals, i.e. integral of \lambda_i until t for each individual
  dt <- object$model$survs[,"stop"] - object$model$survs[,"start"]

  #Easier - no individual summations as they are all summed anyhow afterwards
  intlambda <- tapply(object$model$weights * lambda* dt, BLOCK, sum)

  #Compute cumulative intensities (Ogata (1988): "residual process")
  tau <- cumsum(intlambda)[eventTimesIdx]
  tau
}

