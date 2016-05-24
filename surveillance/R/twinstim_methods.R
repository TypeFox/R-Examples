################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Methods for objects of class "twinstim", specifically:
### vcov, logLik, print, summary, plot, R0, residuals, update, terms, all.equal
###
### Copyright (C) 2009-2016 Sebastian Meyer
### $Revision: 1692 $
### $Date: 2016-04-02 16:24:21 +0200 (Sam, 02. Apr 2016) $
################################################################################

## extract the link function used for the epidemic predictor (default: log-link)
.epilink <- function (x)
{
    link <- attr(x$formula$epidemic, "link")
    if (is.null(link)) "log" else link
}

### don't need a specific coef-method (identical to stats:::coef.default)
## coef.twinstim <- function (object, ...)
## {
##     object$coefficients
## }

## list coefficients by component
coeflist.twinstim <- coeflist.simEpidataCS <- function (x, ...)
{
    coeflist <- coeflist.default(x$coefficients, x$npars)
    ## rename elements and union "nbeta0" and "p" as "endemic"
    coeflist <- c(list(c(coeflist[[1L]], coeflist[[2L]])), coeflist[-(1:2)])
    names(coeflist) <- c("endemic", "epidemic", "siaf", "tiaf")
    coeflist
}

## asymptotic variance-covariance matrix (inverse of expected fisher information)
vcov.twinstim <- function (object, ...)
{
    if (!is.null(object[["fisherinfo"]])) {
        solve(object$fisherinfo)
    } else if (!is.null(object[["fisherinfo.observed"]])) {
        solve(object$fisherinfo.observed)
    } else {
        stop("Fisher information not available; use, e.g., -optimHess()")
    }
}

## Extract log-likelihood of the model (which also enables the use of AIC())
logLik.twinstim <- function (object, ...)
{
    r <- object$loglik
    attr(r, "df") <- length(coef(object))
    attr(r, "nobs") <- nobs(object)
    class(r) <- "logLik"
    r
}

## Also define an extractAIC-method to make step() work
extractAIC.twinstim <- function (fit, scale, k = 2, ...)
{
    loglik <- logLik(fit)
    edf <- attr(loglik, "df")
    penalty <- k * edf
    c(edf = edf, AIC = -2 * c(loglik) + penalty)            
}

## Number of events (excluding the pre-history)
nobs.twinstim <- function (object, ...) length(object$fitted)

## print-method
print.twinstim <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")
    print.default(x$call)
    cat("\nCoefficients:\n")
    print.default(format(coef(x), digits=digits), print.gap = 2, quote = FALSE)
    cat("\nLog-likelihood: ", format(logLik(x), digits=digits), "\n", sep = "")
    if (!isTRUE(x$converged)) {
        cat("\nWARNING: OPTIMIZATION ROUTINE DID NOT CONVERGE!",
            paste0("(",x$converged,")"), "\n")
    }
    cat("\n")
    invisible(x)
}

summary.twinstim <- function (object, test.iaf = FALSE,
    correlation = FALSE, symbolic.cor = FALSE, runtime = FALSE, ...)
{
    ans <- unclass(object)[c("call", "converged", if (runtime) "counts")]
    npars <- object$npars
    nbeta0 <- npars[1]; p <- npars[2]; nbeta <- nbeta0 + p
    q <- npars[3]
    nNotIaf <- nbeta + q
    niafpars <- npars[4] + npars[5]
    est <- coef(object)
    ans$cov <- tryCatch(vcov(object), error = function (e) {
        warning(e)
        matrix(NA_real_, length(est), length(est))
    })
    se <- sqrt(diag(ans$cov))
    zval <- est/se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    coefficients <- cbind(est, se, zval, pval)
    dimnames(coefficients) <- list(names(est),
        c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    ans$coefficients.beta <- coefficients[seq_len(nbeta),,drop=FALSE]
    ans$coefficients.gamma <- structure(
        coefficients[nbeta+seq_len(q),,drop=FALSE],
        link = .epilink(object)
    )
    ans$coefficients.iaf <- coefficients[nNotIaf+seq_len(niafpars),,drop=FALSE]
    if (!test.iaf) {
        ## usually, siaf and tiaf parameters are strictly positive,
        ## or parametrized on the logscale. In this case the usual wald test
        ## with H0: para=0 is invalid or meaningless.
        is.na(ans$coefficients.iaf[,3:4]) <- TRUE
    }
    # estimated parameter correlation
    if (correlation) {
        ans$correlation <- cov2cor(ans$cov)
        ans$symbolic.cor <- symbolic.cor
    }
    ans$loglik <- logLik(object)
    ans$aic <- AIC(object)
    if (runtime) {
        ans$runtime <- object$runtime
    }
    class(ans) <- "summary.twinstim"
    ans
}

## additional methods to make confint.default work for summary.twinstim
vcov.summary.twinstim <- function (object, ...) object$cov
coef.summary.twinstim <- function (object, ...) with(object, {
    coeftab <- rbind(coefficients.beta, coefficients.gamma, coefficients.iaf)
    structure(coeftab[,1], names=rownames(coeftab))
})

## print-method for summary.twinstim
print.summary.twinstim <- function (x,
    digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
    signif.stars = getOption("show.signif.stars"), ...)
{
    nbeta <- nrow(x$coefficients.beta) # = nbeta0 + p
    q <- nrow(x$coefficients.gamma)
    niafpars <- nrow(x$coefficients.iaf)
    cat("\nCall:\n")
    print.default(x$call)
    if (nbeta > 0L) {
        cat("\nCoefficients of the endemic component:\n")
        printCoefmat(x$coefficients.beta, digits = digits,
            signif.stars = signif.stars, signif.legend = (q==0L) && signif.stars, ...)
    } else cat("\nNo coefficients in the endemic component.\n")
    if (q + niafpars > 0L) {
        cat("\nCoefficients of the epidemic component",
            if (attr(x$coefficients.gamma, "link") != "log")
                paste0(" (LINK FUNCTION: ",
                       attr(x$coefficients.gamma, "link"), ")"),
            ":\n", sep = "")
        printCoefmat(rbind(x$coefficients.gamma, x$coefficients.iaf), digits = digits,
            signif.stars = signif.stars, ...)
    } else cat("\nNo epidemic component.\n")
    cat("\nAIC: ", format(x$aic, digits=max(4, digits+1)))
    cat("\nLog-likelihood:", format(x$loglik, digits = digits))
    runtime <- x$runtime
    if (!is.null(runtime)) {
        cat("\nNumber of log-likelihood evaluations:", x$counts[1L])
        cat("\nNumber of score function evaluations:", x$counts[2L])
        cores <- attr(runtime, "cores")
        elapsed <- if (length(runtime) == 1L) { # surveillance < 1.6-0
            runtime
        } else {
            runtime[["elapsed"]]
        }
        cat("\nRuntime",
            if (!is.null(cores) && cores > 1) paste0(" (", cores, " cores)"),
            ": ", format(elapsed, digits = max(4, digits+1)), " seconds",
            sep = "")
    }
    cat("\n")
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
            correl <- format(round(correl, 2), nsmall = 2)
            correl[!lower.tri(correl)] <- ""
            colnames(correl) <- substr(colnames(correl), 1, 5)
            print(correl[-1, -p, drop = FALSE], quote = FALSE)
        }
        }
    }
    if (!isTRUE(x$converged)) {
        cat("\nWARNING: OPTIMIZATION ROUTINE DID NOT CONVERGE!",
            paste0("(",x$converged,")"), "\n")
    }
    cat("\n")
    invisible(x)
}



### 'cat's the summary in LaTeX code

toLatex.summary.twinstim <- function (
    object, digits = max(3, getOption("digits") - 3), eps.Pvalue = 1e-4,
    align = "lrrrr", booktabs = getOption("xtable.booktabs", FALSE),
    withAIC = FALSE, ...)
{
ret <- capture.output({
    cat("\\begin{tabular}{", align, "}\n",
        if (booktabs) "\\toprule" else "\\hline", "\n", sep="")
    cat(" & Estimate & Std. Error & $z$ value & $P(|Z|>|z|)$ \\\\\n",
        if (!booktabs) "\\hline\n", sep="")
    tabh <- object$coefficients.beta
    tabe <- rbind(object$coefficients.gamma, object$coefficients.iaf)
    for (tabname in c("tabh", "tabe")) {
        tab <- get(tabname)
        if (nrow(tab) > 0L) {
            rownames(tab) <- gsub(" ", "", rownames(tab))
            tab_char <- capture.output(
                        printCoefmat(tab, digits=digits, signif.stars=FALSE,
                                     eps.Pvalue = eps.Pvalue, na.print="NA")
                        )[-1]
            ## remove extra space (since used as column sep in read.table)
            tab_char <- sub("< ", "<", tab_char, fixed=TRUE) # small p-values
            ## replace scientific notation by corresponding LaTeX code
            tab_char <- sub("(<?)([0-9]+)e([+-][0-9]+)$",
                            "\\1\\2\\\\cdot{}10^{\\3}",
                            tab_char)
            con <- textConnection(tab_char)
            tab2 <- read.table(con, colClasses="character")
            close(con)
            rownames(tab2) <- paste0("\\texttt{",gsub("_","\\\\_",tab2[,1]),"}")
            tab2 <- tab2[,-1]
            tab2[] <- lapply(tab2, function(x) {
                ifelse(is.na(x), "", paste0("$",x,"$")) # (test.iaf=FALSE)
            })
            cat(if (booktabs) "\\midrule" else "\\hline", "\n")
            print(xtable(tab2), only.contents=TRUE, hline.after=NULL,
                  include.colnames=FALSE, sanitize.text.function=identity)
        }
    }
    if (withAIC) {
        cat(if (booktabs) "\\midrule" else "\\hline", "\n")
        cat("AIC:& $", format(object$aic, digits=max(4, digits+1)), "$ &&&\\\\\n")
        cat("Log-likelihood:& $", format(object$loglik, digits=digits), "$ &&&\\\\\n")
    }
    cat(if (booktabs) "\\bottomrule" else "\\hline", "\n")
    cat("\\end{tabular}\n")
})
class(ret) <- "Latex"
ret
}


## Alternative implementation with exp() of parameters, i.e., rate ratios (RR)
## NOTE: intercepts and iaf parameters are ignored here

xtable.summary.twinstim <- function (x, caption = NULL, label = NULL,
                             align = c("l", "r", "r", "r"), digits = 3,
                             display = c("s", "f", "s", "s"), ...,
                             ci.level = 0.95, ci.fmt = "%4.2f", ci.to = "--",
                             eps.Pvalue = 1e-4)
{
    cis <- confint(x, level=ci.level)
    tabh <- x$coefficients.beta
    tabe <- x$coefficients.gamma
    if (attr(tabe, "link") != "log" && any(rownames(tabe) != "e.(Intercept)"))
        stop("only implemented for the standard log-link models")
    tab <- rbind(tabh, tabe)
    tab <- tab[grep("^([he]\\.\\(Intercept\\)|h.type)", rownames(tab),
                    invert=TRUE),,drop=FALSE]
    expcis <- exp(cis[rownames(tab),,drop=FALSE])
    cifmt <- paste0(ci.fmt, ci.to, ci.fmt)
    rrtab <- data.frame(RR = exp(tab[,1]),
                        CI = sprintf(cifmt, expcis[,1], expcis[,2]),
                        "p-value" = formatPval(tab[,4], eps=eps.Pvalue),
                        check.names = FALSE, stringsAsFactors=FALSE)
    names(rrtab)[2] <- paste0(100*ci.level, "% CI")

    ## append caption etc.
    class(rrtab) <- c("xtable", "data.frame")
    caption(rrtab) <- caption
    label(rrtab) <- label
    align(rrtab) <- align
    digits(rrtab) <- digits
    display(rrtab) <- display

    ## Done
    rrtab
}

xtable.twinstim <- function () {
    cl <- match.call()
    cl[[1]] <- as.name("xtable.summary.twinstim")
    cl$x <- substitute(summary(x))
    eval.parent(cl)                # => xtable.summary.twinstim must be exported
}
formals(xtable.twinstim) <- formals(xtable.summary.twinstim)



### Plot method for twinstim (wrapper for iafplot and intensityplot)

plot.twinstim <- function (x, which, ...)
{
    cl <- match.call()
    which <- match.arg(which, choices =
                       c(eval(formals(intensityplot.twinstim)$which),
                         eval(formals(iafplot)$which)))
    FUN <- if (which %in% eval(formals(intensityplot.twinstim)$which))
        "intensityplot" else "iafplot"
    cl[[1]] <- as.name(FUN)
    if (FUN == "iafplot") names(cl)[names(cl) == "x"] <- "object"
    eval(cl, envir = parent.frame())
}



### Calculates the basic reproduction number R0 for individuals
### with marks given in 'newevents'

R0.twinstim <- function (object, newevents, trimmed = TRUE, newcoef = NULL, ...)
{
    ## check for epidemic component
    npars <- object$npars
    if (npars["q"] == 0L) {
        message("no epidemic component in model, returning 0-vector")
        if (missing(newevents)) return(object$R0) else {
            return(structure(rep.int(0, nrow(newevents)),
                             names = rownames(newevents)))
        }
    }
    ## update object for use of new parameters
    if (!is.null(newcoef)) {
        object <- update(object, optim.args = list(par=newcoef, fixed=TRUE),
                         cumCIF = FALSE, cores = 1L, verbose = FALSE)
    }
    ## extract model information
    t0 <- object$timeRange[1L]
    T <- object$timeRange[2L]
    typeNames <- rownames(object$qmatrix)
    nTypes <- length(typeNames)
    types <- seq_len(nTypes)
    form <- formula(object)
    siaf <- form$siaf
    tiaf <- form$tiaf
    coefs <- coef(object)
    tiafpars <- coefs[sum(npars[1:4]) + seq_len(npars["ntiafpars"])]
    siafpars <- coefs[sum(npars[1:3]) + seq_len(npars["nsiafpars"])]
    
    if (missing(newevents)) {
        ## if no newevents are supplied, use original events
        if (trimmed) {                  # already calculated by 'twinstim'
            return(object$R0)
        } else {    # untrimmed version (spatio-temporal integral over R+ x R^2)
            ## extract relevant data from model environment
            if (is.null(modelenv <- environment(object))) {
                stop("need model environment for untrimmed R0 of fitted events\n",
                     " -- re-fit or update() with 'model=TRUE'")
            }
            eventTypes <- modelenv$eventTypes
            eps.t <- modelenv$eps.t
            eps.s <- modelenv$eps.s
            gammapred <- modelenv$gammapred
            names(gammapred) <- names(object$R0) # for names of the result
        }
    } else {
        ## use newevents
        stopifnot(is.data.frame(newevents))
        if (!"time" %in% names(newevents)) {
            stop("missing event \"time\" column in 'newevents'")
        }
        if (any(!c("eps.s", "eps.t") %in% names(newevents))) {
            stop("missing \"eps.s\" or \"eps.t\" columns in 'newevents'")
        }
        stopifnot(is.factor(newevents[["type"]]))
        
        ## subset newevents to timeRange
        .N <- nrow(newevents)
        newevents <- subset(newevents, time + eps.t > t0 & time <= T)
        if (nrow(newevents) < .N) {
            message("subsetted 'newevents' to only include events infectious ",
                    "during 'object$timeRange'")
        }

        ## extract columns
        newevents$type <- factor(newevents[["type"]], levels = typeNames)
        eventTimes <- newevents[["time"]]
        eps.t <- newevents[["eps.t"]]
        eps.s <- newevents[["eps.s"]]
        
        ## calculate gammapred for newevents
        epidemic <- terms(form$epidemic, data = newevents, keep.order = TRUE)
        mfe <- model.frame(epidemic, data = newevents,
                           na.action = na.pass, drop.unused.levels = FALSE)
        mme <- model.matrix(epidemic, mfe)
        gamma <- coefs[sum(npars[1:2]) + seq_len(npars["q"])]
        if (ncol(mme) != length(gamma)) {
            stop("epidemic model matrix has the wrong number of columns ",
                 "(check the variable types in 'newevents' (factors, etc))")
        }
        gammapred <- drop(mme %*% gamma)  # identity link
        if (.epilink(object) == "log")
            gammapred <- exp(gammapred)
        names(gammapred) <- rownames(newevents)

        ## now, convert types of newevents to integer codes
        eventTypes <- as.integer(newevents$type)
    }

    ## qSum
    qSumTypes <- rowSums(object$qmatrix)
    qSum <- unname(qSumTypes[eventTypes])


    ## calculate remaining factors of the R0 formula, i.e. siafInt and tiafInt
    
    if (trimmed) {                      # trimmed R0 for newevents
        
        ## integral of g over the observed infectious periods
        .tiafInt <- .tiafIntFUN()
        gIntUpper <- pmin(T - eventTimes, eps.t)
        gIntLower <- pmax(0, t0 - eventTimes)
        tiafInt <- .tiafInt(tiafpars, from=gIntLower, to=gIntUpper,
                            type=eventTypes, G=tiaf$G)
        ## integral of f over the influenceRegion
        bdist <- newevents[[".bdist"]]
        influenceRegion <- newevents[[".influenceRegion"]]
        if (is.null(influenceRegion)) {
            stop("missing \".influenceRegion\" component in 'newevents'")
        }
        noCircularIR <- if (is.null(bdist)) FALSE else all(eps.s > bdist)
        if (attr(siaf, "constant")) {
            iRareas <- sapply(influenceRegion, area.owin)
            ## will be used by .siafInt()
        } else if (! (is.null(siaf$Fcircle) ||
               (is.null(siaf$effRange) && noCircularIR))) {
            if (is.null(bdist)) {
                stop("missing \".bdist\" component in 'newevents'")
            }
        }
        .siafInt <- .siafIntFUN(siaf, noCircularIR=noCircularIR)
        .siafInt.args <- c(alist(siafpars), object$control.siaf$F)
        siafInt <- do.call(".siafInt", .siafInt.args)
        
    } else {                     # untrimmed R0 for original events or newevents

        ## integrals of interaction functions for all combinations of type and
        ## eps.s/eps.t in newevents
        typeTcombis <- expand.grid(type=types, eps.t=unique(eps.t),
                                   KEEP.OUT.ATTRS=FALSE)
        typeTcombis$gInt <-
            with(typeTcombis, tiaf$G(eps.t, tiafpars, type)) -
                tiaf$G(rep.int(0,nTypes), tiafpars, types)[typeTcombis$type]

        Fcircle <- getFcircle(siaf, object$control.siaf$F)
        typeScombis <- expand.grid(type=types, eps.s=unique(eps.s),
                                   KEEP.OUT.ATTRS=FALSE)
        typeScombis$fInt <- apply(typeScombis, MARGIN=1, FUN=function (type_eps.s) {
            type <- type_eps.s[1L]
            eps.s <- type_eps.s[2L]
            Fcircle(eps.s, siafpars, type)
        })
        
        ## match combinations to rows of original events or 'newevents'
        eventscombiidxS <- match(paste(eventTypes,eps.s,sep="."),
                                 with(typeScombis,paste(type,eps.s,sep=".")))
        eventscombiidxT <- match(paste(eventTypes,eps.t,sep="."),
                                 with(typeTcombis,paste(type,eps.t,sep=".")))

        siafInt <- typeScombis$fInt[eventscombiidxS]
        tiafInt <- typeTcombis$gInt[eventscombiidxT]

        if (any(is.infinite(eps.t) & !is.finite(tiafInt),
                is.infinite(eps.s) & !is.finite(siafInt))) {
            message("infinite interaction ranges yield non-finite R0 values ",
                    "because 'trimmed = FALSE'")
        }
        
    }

    ## return R0 values
    R0s <- qSum * gammapred * siafInt * tiafInt
    R0s
}

## calculate simple R0 (over circular domain, without epidemic covariates,
## for type-invariant siaf/tiaf)
simpleR0 <- function (object, eta = coef(object)[["e.(Intercept)"]],
                      eps.s = NULL, eps.t = NULL, newcoef = NULL)
{
    stopifnot(inherits(object, c("twinstim", "simEpidataCS")))
    if (object$npars[["q"]] == 0L)
        return(0)
    if (any(rowSums(object$qmatrix) != 1))
        warning("'simpleR0' is not correct for type-specific epidemic models")

    if (!is.null(newcoef)) { # use alternative coefficients
        object$coefficients <- newcoef
    }
    coeflist <- coeflist(object)
    siaf <- object$formula$siaf
    tiaf <- object$formula$tiaf

    ## default radii of interaction
    if (is.null(eps.s)) {
        eps.s <- attr(siaf, "eps")
        if (length(eps.s) > 1L) stop("found non-unique 'eps.s'; please set one")
    } else stopifnot(isScalar(eps.s))
    if (is.null(eps.t)) {
        eps.t <- attr(tiaf, "eps")
        if (length(eps.t) > 1L) stop("found non-unique 'eps.t'; please set one")
    } else stopifnot(isScalar(eps.t))

    ## integral of siaf over a disc of radius eps.s
    Fcircle <- getFcircle(siaf, object$control.siaf$F)
    siafInt <- Fcircle(eps.s, coeflist$siaf)

    ## integral of tiaf over a period of length eps.t
    tiafInt <- tiaf$G(eps.t, coeflist$tiaf) - tiaf$G(0, coeflist$tiaf)

    ## calculate basic R0
    (if (.epilink(object) == "log") exp(eta) else eta) * siafInt * tiafInt
}            



### Extract the "residual process" (cf. Ogata, 1988) of a twinstim, i.e. the
### fitted cumulative intensity of the ground process at the event times.
### "generalized residuals similar to those discussed in Cox and Snell (1968)"

residuals.twinstim <- function (object, ...)
{
  res <- object$tau
  if (is.null(res)) {
      if (is.null(modelenv <- environment(object))) {
          stop("residuals not available; re-fit the model with 'cumCIF = TRUE'")
      } else {
          message("'", substitute(object), "' was fit with disabled 'cumCIF'",
                  " -> calculate it now ...")
          res <- with(modelenv, LambdagEvents(cumCIF.pb = interactive()))
          try({
              objname <- deparse(substitute(object))
              object$tau <- res
              assign(objname, object, envir = parent.frame())
              message("Note: added the 'tau' component to object '", objname,
                      "' for future use.")
          }, silent = TRUE)
      }
  }
  return(res)
}



######################################################################
# Function to compute estimated and profile likelihood based
# confidence intervals. Heavy computations might be necessary!
#
#Params:
# fitted - output from a fit with twinstim
# profile - list with 4D vector as entries - format:
#               c(index, lower, upper, grid size)
#           where index is the index in the coef vector
#                 lower and upper are the parameter limits (can be NA)
#                 grid size is the grid size of the equally spaced grid
#                 between lower and upper (can be 0)
# alpha - (1-alpha)% profile likelihood CIs are computed.
#         If alpha <= 0 then no CIs are computed
# control - control object to use for optim in the profile loglik computations
#
# Returns:
#  list with profile loglikelihood evaluations on the grid
#  and highest likelihood and wald confidence intervals
######################################################################

profile.twinstim <- function (fitted, profile, alpha = 0.05,
    control = list(fnscale = -1, maxit = 100, trace = 1),
    do.ltildeprofile=FALSE, ...)
{
  warning("the profile likelihood implementation is experimental")
  ## the implementation below is not well tested, simply uses optim (ignoring
  ## optimizer settings from the original fit), and does not store the complete
  ## set of coefficients
  
  ## Check that input is ok
  profile <- as.list(profile)
  if (length(profile) == 0L) {
    stop("nothing to do")
  }
  lapply(profile, function(one) {
    if (length(one) != 4L) {
      stop("each profile entry has to be of form ",
           "'c(index, lower, upper, grid size)'")
    }})
  if (is.null(fitted[["functions"]])) {
    stop("'fitted' must contain the component 'functions' -- fit using the option model=TRUE")
  }

  ## Control of the optim procedure
  if (is.null(control[["fnscale",exact=TRUE]])) { control$fnscale <- -1 }
  if (is.null(control[["maxit",exact=TRUE]])) { control$maxit <- 100 }
  if (is.null(control[["trace",exact=TRUE]])) { control$trace <- 1 }


  ## Estimated normalized likelihood function
  ltildeestim <- function(thetai,i) {
    theta <- theta.ml
    theta[i] <- thetai
    fitted$functions$ll(theta) - loglik.theta.ml
  }

  ## Profile normalized likelihood function
  ltildeprofile <- function(thetai,i)
  {
    #cat("Investigating theta[",i,"] = ",thetai,"\n")

    emptyTheta <- rep(0, length(theta.ml))

    # Likelihood l(theta_{-i}) = l(theta_i, theta_i)
    ltildethetaminusi <- function(thetaminusi) {
      theta <- emptyTheta
      theta[-i] <- thetaminusi
      theta[i] <- thetai
      #cat("Investigating theta = ",theta,"\n")
      res <- fitted$functions$ll(theta) - loglik.theta.ml
      #cat("Current ltildethetaminusi value: ",res,"\n")
      return(res)
    }
    # Score function of all params except thetaminusi
    stildethetaminusi <- function(thetaminusi) {
      theta <- emptyTheta
      theta[-i] <- thetaminusi
      theta[i] <- thetai
      res <- fitted$functions$sc(theta)[-i]
      #cat("Current stildethetaminusi value: ",res,"\n")
      return(res)
    }

    # Call optim -- currently not adapted to arguments of control arguments
    # used in the fit
    resOthers <- tryCatch(
            optim(par=theta.ml[-i], fn = ltildethetaminusi, gr = stildethetaminusi,
                  method = "BFGS", control = control),
            error = function(e) list(value=NA))
    resOthers$value
  }



  ## Initialize
  theta.ml <- coef(fitted)
  loglik.theta.ml <- c(logLik(fitted))
  se <- sqrt(diag(vcov(fitted)))
  resProfile <- list()


  ## Perform profile computations for all requested parameters
  cat("Evaluating the profile logliks on a grid...\n")
  for (i in 1:length(profile))
    {
    cat("i= ",i,"/",length(profile),"\n")
    #Index of the parameter in the theta vector
    idx <- profile[[i]][1]
    #If no borders are given use those from wald intervals (unconstrained)
    if (is.na(profile[[i]][2])) profile[[i]][2] <- theta.ml[idx] - 3*se[idx]
    if (is.na(profile[[i]][3])) profile[[i]][3] <- theta.ml[idx] + 3*se[idx]
    #Evaluate profile loglik on a grid (if requested)
    if (profile[[i]][4] > 0) {
      thetai.grid <- seq(profile[[i]][2],profile[[i]][3],length=profile[[i]][4])
      resProfile[[i]] <- matrix(NA, nrow = length(thetai.grid), ncol = 4L,
        dimnames = list(NULL, c("grid","profile","estimated","wald")))

      #Loop over all gridpoints
      for (j in 1:length(thetai.grid)) {
        cat("\tj= ",j,"/",length(thetai.grid),"\n")
        resProfile[[i]][j,] <- c(thetai.grid[j],
           #Do we need to compute ltildeprofile (can be quite time consuming)
           if (do.ltildeprofile) ltildeprofile(thetai.grid[j],idx) else NA_real_,
           ltildeestim(thetai.grid[j],idx),
           - 1/2*(1/se[idx]^2)*(thetai.grid[j] - theta.ml[idx])^2)
      }
    }
  }
  names(resProfile) <- names(theta.ml)[sapply(profile, function(x) x[1L])]

  ###############################
  ## Profile likelihood intervals
  ###############################
  # Not done, yet
  ciProfile <- NULL

  ####Done, return
  return(list(lp=resProfile, ci.hl=ciProfile, profileObj=profile))
}



### update-method for the twinstim-class
## stats::update.default would also work but is not adapted to the specific
## structure of twinstim: optim.args (use modifyList), two formulae, model, ...
## However, this specific method is inspired by and copies small parts of the
## update.default method from the stats package developed by The R Core Team

update.twinstim <- function (object, endemic, epidemic,
                             control.siaf, optim.args, model,
                             ..., use.estimates = TRUE, evaluate = TRUE)
{
    call <- object$call
    thiscall <- match.call(expand.dots=FALSE)
    extras <- thiscall$...
    
    if (!missing(model)) {
        call$model <- model
        ## Special case: update model component ONLY
        if (evaluate &&
            all(names(thiscall)[-1] %in% c("object", "model", "evaluate"))) {
            return(.update.twinstim.model(object, model))
        }
    }

    ## Why we no longer use call$endemic but update object$formula$endemic:
    ## call$endemic would be an unevaluated expression eventually receiving the
    ## parent.frame() as environment, cp.: 
    ##(function(e) {ecall <- match.call()$e; eval(call("environment", ecall))})(~1+start)
    ## This could cause large files if the fitted model is saved.
    ## Furthermore, call$endemic could refer to some object containing
    ## the formula, which is no longer visible.    
    call$endemic <- if (missing(endemic)) object$formula$endemic else
        update.formula(object$formula$endemic, endemic)
    call$epidemic <- if (missing(epidemic)) object$formula$epidemic else
        update.formula(object$formula$epidemic, epidemic)
    ## Note: update.formula uses terms.formula(...,simplify=TRUE), but
    ##       the principle order of terms is retained. Offsets will be moved to 
    ##       the end and a missing intercept will be denoted by a final -1.
    
    if (!missing(control.siaf)) {
        if (is.null(control.siaf)) {
            call$control.siaf <- NULL  # remove from call, i.e., use defaults
        } else {
            call$control.siaf <- object$control.siaf # =NULL if constantsiaf
            call$control.siaf[names(control.siaf)] <- control.siaf
        }
    }
    
    call["optim.args"] <- if (missing(optim.args)) object["optim.args"] else {
        list( # use list() to enable optim.args=NULL
             if (is.list(optim.args)) {
                 modifyList(object$optim.args, optim.args)
             } else optim.args           # = NULL
             )
    }
    ## Set initial values (will be appropriately subsetted and/or extended with
    ## zeroes inside twinstim())
    call$start <- if (missing(optim.args) ||
                      (!is.null(optim.args) && !"par" %in% names(optim.args))) {
        ## old optim.args$par probably doesn't match updated model,
        ## thus we set it as "start"-argument
        call$optim.args$par <- NULL
        if (use.estimates) coef(object) else object$optim.args$par
    } else NULL
    if ("start" %in% names(extras)) {
        newstart <- check_twinstim_start(eval.parent(extras$start))
        call$start[names(newstart)] <- newstart
        extras$start <- NULL
    }
    ## CAVE: the remainder is copied from stats::update.default (as at R-2.15.0)
    if(length(extras)) {
	existing <- !is.na(match(names(extras), names(call)))
	## do these individually to allow NULL to remove entries.
	for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
	if(any(!existing)) {
	    call <- c(as.list(call), extras[!existing])
	    call <- as.call(call)
	}
    }
    if(evaluate) eval(call, parent.frame())
    else call
}

.update.twinstim.model <- function (object, model)
{
    call <- object$call
    call$model <- model
    if (model) { # add model environment
        call$start <- coef(object)
        call$optim.args$fixed <- TRUE
        call$cumCIF <- FALSE
        call$verbose <- FALSE
        ## evaluate in the environment calling update.twinstim()
        message("Setting up the model environment ...")
        objectWithModel <- eval(call, parent.frame(2L))
        ## add the model "functions" and environment
        object$functions <- objectWithModel$functions
        environment(object) <- environment(objectWithModel)
    } else { # remove model environment
        object["functions"] <- list(NULL)
        environment(object) <- NULL
    }
    object$call$model <- model
    object
}

## a terms-method is required for stepComponent()
terms.twinstim <- function (x, component=c("endemic", "epidemic"), ...)
{
    component <- match.arg(component)
    terms.formula(x$formula[[component]], keep.order=TRUE)
}

## compare two twinstim fits ignoring at least the "runtime" and the "call"
## just like all.equal.hhh4()
all.equal.twinstim <- function (target, current, ..., ignore = NULL)
{
    ignore <- unique.default(c(ignore, "runtime", "call"))
    target[ignore] <- current[ignore] <- list(NULL)
    NextMethod("all.equal")
}
