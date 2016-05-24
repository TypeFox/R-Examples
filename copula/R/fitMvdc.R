## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

setClass("fitMvdc", representation(mvdc = "mvdc"),
	 contains="fittedMV" #-> ./Classes.R
	 ## FIXME , validity = function(object) TRUE
	 )

paramNmsMvdc <- function(mv) c(margpnames(mv), mv@copula@param.names)
## not exported, but used ..
setMethod("paramNames", "fitMvdc", function(x) paramNmsMvdc(x@mvdc))

print.fitMvdc <- function(x, digits = max(3, getOption("digits") - 3),
                          signif.stars = getOption("show.signif.stars"), ...)
{
    foo <- summary.fitMvdc(x)
    cat("The Maximum Likelihood estimation is based on", x@nsample, "observations.\n")
    p <- x@mvdc@copula@dimension
    marNpar <- vapply(x@mvdc@paramMargins, length, 1L)
    idx2 <- cumsum(marNpar)
    idx1 <- idx2 - marNpar + 1
    margid <- x@mvdc@marginsIdentical
    if (sum(marNpar) > 0) { ## sometimes there is no marginal params
	if(margid){
	    cat("Identical margins:\n")
	    printCoefmat(foo$coefficients[idx1[1]:idx2[1], , drop=FALSE],
			 digits = digits, signif.stars = signif.stars,
			 signif.legend=FALSE, ## in any case
			 na.print = "NA", ...)
	}
	else {
	    for (i in 1:p) {
		cat("Margin", i, ":\n")
		printCoefmat(foo$coefficients[idx1[i]:idx2[i], 1:2, drop=FALSE],
			     digits = digits, signif.stars = signif.stars,
			     signif.legend=FALSE, ## in any case
			     na.print = "NA", ...)

	    }
	}
    }
    cat("Copula:\n")
    copParIdx <- seq_along(x@mvdc@copula@parameters)
    printCoefmat(foo$coefficients[copParIdx +
				  if(margid) marNpar[1] else sum(marNpar),
				  if(margid) 1:4 else 1:2, drop=FALSE],
		 digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    cat("The maximized loglikelihood is", foo$loglik, "\n")
    if (!is.na(foo$convergence)) {
	if(foo$convergence)
	    cat("Convergence problems: code is", foo$convergence, "see ?optim.\n")
	else cat("Optimization converged\n")
    }
    if(!is.null(cnts <- x@fitting.stats$counts) && !all(is.na(cnts))) {
	cat("Number of loglikelihood evaluations:\n"); print(cnts, ...)
    }
    invisible(x)
}

summary.fitMvdc <- function(object, ...) {
  estimate <- object@estimate
  se <- sqrt(diag(object@var.est))
  ## Hmm, FIXME?  'theta=0' is often *not* the interesting null-hypothesis
  ## neither for the margins, nor the copula
  zvalue <- estimate / se
  pvalue <- 2* pnorm(abs(zvalue), lower.tail=FALSE)
  ##ans <- object[c("loglik", "convergence")]
  coef <- cbind(Estimate = estimate, "Std. Error" = se,
                "z value" = zvalue, "Pr(>|z|)" = pvalue)
  rownames(coef) <- paramNmsMvdc(object@mvdc)
  structure(class = "summary.fitMvdc",
	    list(loglik = object@loglik,
		 convergence = object@fitting.stats[["convergence"]],
		 coefficients = coef))
}

setMethod("show", signature("fitMvdc"), function(object) print.fitMvdc(object))


################################################################################

setMvdcPar <- function(mvdc, param) {
    marNpar <- vapply(mvdc@paramMargins, length, 1L)
    idx2 <- cumsum(marNpar)
    idx1 <- idx2 - marNpar + 1
    margid <- mvdc@marginsIdentical
    p <- mvdc@copula@dimension

    for (i in 1:p) {
        if (marNpar[i] > 0) {
            ## parnames <- mvdc@paramMargins[[i]]
            k <- if(margid) 1 else i
            par <- param[idx1[k]: idx2[k]]
            ## names(par) <- parnames
            ## mvdc@paramMargins[i] <- as.list(par)
            for (j in 1:marNpar[i]) mvdc@paramMargins[[i]][j] <- par[j]
        }
    }
    mvdc@copula@parameters <-
        if (idx2[p] == 0)               # no marginal parameters
            param else param[- (if(margid) 1:idx2[1] else 1:rev(idx2)[1])]
    mvdc
}

loglikMvdc <- function(param, x, mvdc, hideWarnings) {
  if(!missing(hideWarnings))
      warning("'hideWarnings' is deprecated and has no effect here anymore")
  ##       use suppressMessages() otherwise {and live without  "fitMessages"}

  mvdc <- setMvdcPar(mvdc, param)

  loglik <- tryCatch(sum(log(dMvdc(x, mvdc))), error = function(e) e)

  if(is(loglik, "error")) {
      warning("error in loglik computation: ", loglik$message)
      (-Inf)# was NaN
  }
  else loglik
}

fitMvdc <- function(data, mvdc, start,
                    optim.control=list(), method="BFGS",
                    lower = -Inf, upper = Inf,
                    estimate.variance = fit$convergence == 0, hideWarnings=TRUE)
{
    copula <- mvdc@copula
    if (copula@dimension != ncol(data))
        stop("The dimensions of the data and copula do not match.")
    marNpar <- vapply(mvdc@paramMargins, length, 1L)
    margid <- mvdc@marginsIdentical
    q <- length(start)
    if(q != length(copula@parameters) + (if(margid) marNpar[1] else sum(marNpar)))
	stop("The lengths of 'start' and mvdc parameters do not match.")
    mvdCheckM(mvdc@margins, "p")
    control <- c(optim.control, fnscale=-1)
    control <- control[ !vapply(control, is.null, NA)]

    ## messageOut may be used for debugging
    if (hideWarnings) {
	messageOut <- textConnection("fitMessages", open="w", local=TRUE)
	sink(messageOut); sink(messageOut, type="message")
	oop <- options(warn = -1) ## ignore warnings; can be undesirable!
	on.exit({ options(oop); sink(type="message"); sink(); close(messageOut)})
    }

    fit <- optim(start, loglikMvdc,
		 ## loglikMvdc args :
		 mvdc=mvdc, x=data,
		 ## optim args:
		 method=method, control=control, lower=lower, upper=upper)

    if (hideWarnings) {
	options(oop); sink(type="message"); sink()
	on.exit()
        close(messageOut)
    }

    if (fit$convergence > 0)
	warning("possible convergence problem: optim gave code=", fit$convergence)
    loglik <- fit$val
    param <- fit$par

    varNA <- matrix(NA_real_, q, q)
    var.est <- if (estimate.variance) {
	fit.last <- optim(param, loglikMvdc,
			  ## loglikMvdc args :
			  mvdc=mvdc, x=data,
			  ## optim args:
			  method = method, ## one final step, computing Hessian :
			  control=c(control, maxit = 1), hessian=TRUE)
        fit$counts <- fit$counts + fit.last$counts
	vcov <- tryCatch(solve(-fit.last$hessian), error = function(e) e)
	if(is(vcov, "error")) {
	    warning("Hessian matrix not invertible: ", vcov$message)
	    varNA
	} else vcov ## ok
    } else varNA

    new("fitMvdc",
	estimate = param,
	var.est = var.est,
	loglik = loglik,
	method = method,
	fitting.stats = c(fit[c("convergence", "counts", "message")], control),
	nsample = nrow(data),
	## this contains 'copula':
	mvdc = setMvdcPar(mvdc, param))
}
