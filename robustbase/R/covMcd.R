### This is originally from the R package
####
####  rrcov : Scalable Robust Estimators with High Breakdown Point
####
#### by Valentin Todorov

##  I would like to thank Peter Rousseeuw and Katrien van Driessen for
##  providing the initial code of this function.

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## http://www.r-project.org/Licenses/

## No longer hidden in namespace :
## easier to explain when user-available & documented if
h.alpha.n <- function(alpha, n, p) {
    ## Compute h(alpha) := size of subsample, given alpha, (n,p)
    ## Same function for covMcd() and ltsReg()
    n2 <- (n+p+1) %/% 2
    floor(2 * n2 - n + 2 * (n - n2) * alpha)
}

## MM: the way it's set up, *must* be kept in sync with rrcov.control()'s
## defaults --> ./rrcov.control.R  :
covMcd <- function(x,
           cor = FALSE,
	   raw.only = FALSE,
           alpha = control$ alpha,
           nsamp = control$ nsamp,
           nmini = control$ nmini, kmini = control$ kmini,
           scalefn=control$scalefn, maxcsteps=control$maxcsteps,
           initHsets = NULL, save.hsets = FALSE, names = TRUE,
           seed  = control$ seed,
           tolSolve = control$ tolSolve, # had 1e-10 hardwired {now 1e-14 default}
           trace = control$ trace,
           use.correction = control$ use.correction,
           wgtFUN = control$ wgtFUN,
           control = rrcov.control())
{
    logdet.Lrg <- 50 ## <-- FIXME add to  rrcov.control() and then use that
    ##   Analyze and validate the input parameters ...
    if(length(seed) > 0) {
	if(length(seed) < 3 || seed[1L] < 100)
	    stop("invalid 'seed'. Must be compatible with .Random.seed !")
        if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
            seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
        }
        assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    ## For back compatibility, as some new args did not exist pre 2013-04,
    ## and callers of covMcd() may use a "too small"  'control' list:
    defCtrl <- if(missing(control)) control else rrcov.control()
    if(missing(wgtFUN)) getDefCtrl("wgtFUN", defCtrl)
    if(is.null (nmini)) getDefCtrl("nmini", defCtrl)

    ##   vt::03.02.2006 - added options "best" and "exact" for nsamp
    ##   nsamp will be further analized in the wrapper .fastmcd()
    if(is.numeric(nsamp) && nsamp <= 0)
        stop("Invalid number of trials nsamp = ",nsamp, "!")

    if(is.data.frame(x))
	x <- data.matrix(x, rownames.force=FALSE)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
                    dimnames = if(names) list(names(x), deparse(substitute(x))))

    if(!names) dimnames(x) <- NULL # (speedup)
    ## drop all rows with missing values (!!) :
    ok <- is.finite(x %*% rep.int(1, ncol(x)))
    x <- x[ok, , drop = FALSE]
    if(!length(dx <- dim(x)))
        stop("All observations have missing values!")
    n <- dx[1]; p <- dx[2]
    if(names) dimn <- dimnames(x)
    ## h(alpha) , the size of the subsamples
    h <- h.alpha.n(alpha, n, p)
    if(n <= p + 1) # ==> floor((n+p+1)/2) > n - 1  -- not Ok
        stop(if (n <= p) # absolute barrier!
             "n <= p -- you can't be serious!"
        else "n == p+1  is too small sample size for MCD")
    ## else
    if(n < 2 * p) { ## p+1 < n < 2p
        warning("n < 2 * p, i.e., possibly too small sample size")
        ## was stop("Need at least 2*(number of variables) observations ")
    }
    ##     jmin <- (n + p + 1) %/% 2
    ##     if(alpha < 1/2) ## FIXME? shouldn't we rather test   'alpha < jmin/n' ?
    ##  stop("The MCD must cover at least", jmin, "observations")
    ## MM: I think this should be sufficient;
    ##     we should even omit the (n < 2p) warning
    if(h > n)
        stop("Sample size n  <  h(alpha; n,p) := size of \"good\" subsample")
    else if(2*h < n)
	warning("subsample size	 h < n/2  may be too small")

    if(is.character(wgtFUN)) {
	if(is.function(mkWfun <- .wgtFUN.covMcd[[wgtFUN]]))
            wgtFUN <- mkWfun(p=p, n=n, control)
    }
    if(!is.function(wgtFUN))
	stop(gettextf("'wgtFUN' must be a function or one of the strings %s.",
		      pasteK(paste0('"',names(.wgtFUN.covMcd),'"'))), domain=NA)

    ## vt::03.02.2006 - raw.cnp2 and cnp2 are vectors of size 2 and  will
    ##   contain the correction factors (concistency and finite sample)
    ##   for the raw and reweighted estimates respectively. Set them
    ##   initially to 1.  If use.correction is false (not the default),
    ##   the finite sample correction factor will not be used
    ##   (neither for the raw estimates nor for the reweighted ones)
    raw.cnp2 <- cnp2 <- c(1,1)

    ans <- list(call = match.call(), nsamp = nsamp,
                method = sprintf("MCD(alpha=%g ==> h=%d)", alpha, h))

    if(h == n) { ## <==> alpha ~= 1 : Just compute the classical estimates --------
        mcd <- cov(x) #MM: was  cov.wt(x)$cov
        loc <- as.vector(colMeans(x))
        obj <- determinant(mcd, logarithm = TRUE)$modulus[1]
        if ( -obj/p > logdet.Lrg ) {
            ans$cov <- mcd
	    if(names) dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
            if (cor)
                ans$cor <- cov2cor(ans$cov)
            ans$center <- loc
            if(names && length(dimn[[2]]))
                names(ans$center) <- dimn[[2]]
            ans$n.obs <- n
            ans$singularity <- list(kind = "classical")
            weights <- 1
        }
        else {
            mah <- mahalanobis(x, loc, mcd, tol = tolSolve)
            ## VT:: 01.09.2004 - bug in alpha=1
            weights <- wgtFUN(mah) # 0/1
            sum.w <- sum(weights)
            ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
            ## cov.wt() -> list("cov", "center", "n.obs", ["wt", "cor"])
            ## Consistency factor for reweighted MCD -- ok for default wgtFUN only: FIXME
            if(sum.w != n) {
                cnp2[1] <- .MCDcons(p, sum.w/n)
                ans$cov <- ans$cov * cnp2[1]
            }
            obj <- determinant(mcd, logarithm = TRUE)$modulus[1]
            if( -obj/p > logdet.Lrg ) {
                ans$singularity <- list(kind = "reweighted.MCD")
            }
            else {
                mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
                weights <- wgtFUN(mah) # 0/1
            }
        }

        ans$alpha <- alpha
        ans$quan <- h
        ans$raw.cov <- mcd
        ans$raw.center <- loc
        if(names && !is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$crit <- obj # was exp(obj); but log-scale is "robust" against under/overflow
        ans$method <- paste(ans$method,
                            "\nalpha = 1: The minimum covariance determinant estimates based on",
                            n, "observations \nare equal to the classical estimates.")
        ans$mcd.wt <- rep.int(NA, length(ok))
        ans$mcd.wt[ok] <- weights
        if(names && length(dimn[[1]]))
            names(ans$mcd.wt) <- dimn[[1]]
        ans$wt <- NULL
        ans$X <- x
        if(names) {
            if(length(dimn[[1]]))
                dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
            else
                dimnames(ans$X) <- list(seq(along = ok)[ok], NULL)
        }
        if(trace)
            cat(ans$method, "\n")
        ans$raw.cnp2 <- raw.cnp2
        ans$cnp2 <- cnp2
        class(ans) <- "mcd"
        return(ans)
    } ## end { alpha = 1   <==>   h = n }

    mcd <- if(nsamp == "deterministic") {
	ans$method <- paste("Deterministic", ans$method)
	.detmcd (x, h, hsets.init = initHsets,
		 save.hsets=save.hsets, # full.h=full.h,
		 scalefn=scalefn, maxcsteps=maxcsteps, trace=as.integer(trace),
		 names=names)
    } else {
	ans$method <- paste0("Fast ", ans$method, "; nsamp = ", nsamp,
			     "; (n,k)mini = (", nmini,",",kmini,")")
	.fastmcd(x, h, nsamp, nmini, kmini, trace=as.integer(trace))
    }

    ## Compute the consistency correction factor for the raw MCD
    ##  (see calfa in Croux and Haesbroeck)
    calpha <- .MCDcons(p, h/n)    ## VT::19.3.2007
    correct <- if(use.correction) .MCDcnp2(p, n, alpha) else 1.
    raw.cnp2 <- c(calpha, correct)

    if(p == 1) {
        ## ==> Compute univariate location and scale estimates
	ans$method <- paste("Univariate", ans$method)
        scale <- sqrt(calpha * correct) * as.double(mcd$initcovariance)
        center <- as.double(mcd$initmean)
        if(abs(scale - 0) < 1e-07) {
            ans$singularity <- list(kind = "identicalObs", q = h)
            ans$raw.cov <- ans$cov <- matrix(0)
            ans$raw.center <- ans$center <- center
            ans$n.obs <- n
            ans$alpha <- alpha
            ans$quan <- h
            if(names && !is.null(nms <- dimn[[2]][1])) {
                names(ans$raw.center) <- names(ans$center) <- nms
                dimnames(ans$raw.cov) <- dimnames(ans$cov) <- list(nms,nms)
            }
            ans$crit <- -Inf # = log(0)
            weights <- as.numeric(abs(x - center) < 1e-07) # 0 / 1
        } ## end { scale ~= 0 }
        else {
            ## Compute the weights for the raw MCD in case p=1
            weights <- wgtFUN(((x - center)/scale)^2) # 0/1
            sum.w <- sum(weights)
            ans <- c(ans, cov.wt(x, wt = weights, cor=cor))

	    if(sum.w != n) {
		cdelta.rew <- .MCDcons(p, sum.w/n) ## VT::19.3.2007
		correct.rew <- if(use.correction) .MCDcnp2.rew(p, n, alpha) else 1.
		cnp2 <- c(cdelta.rew, correct.rew)
		ans$cov <- cdelta.rew * correct.rew * ans$cov
	    }
            ans$alpha <- alpha
            ans$quan <- h
            ans$raw.cov <- as.matrix(scale^2)
            ans$raw.center <- as.vector(center)
            if(names && !is.null(nms <- dimn[[2]][1])) {
                dimnames(ans$raw.cov) <- list(nms,nms)
                names(ans$raw.center) <- nms
            }
	    ans$crit <- ## log(det) =
		log(sum(sort((x - as.double(mcd$initmean))^2, partial = h)[1:h])/max(1,h-1))
            center <- ans$center
            scale <- as.vector(sqrt(ans$cov))
            weights <- wgtFUN(((x - center)/scale)^2)
        } ## end{ scale > 0 }
    } ## end p=1

    else { ## p >= 2 : ---------------------------------------------------------

      ## Apply correction factor to the raw estimates
      ## and use them to compute weights
      mcd$initcovariance <- matrix(calpha * correct * mcd$initcovariance, p,p)
      if(raw.only || mcd$exactfit != 0) {
        ## If not all observations are in general position, i.e. more than
        ## h observations lie on a hyperplane, the program still yields
        ## the MCD location and scatter matrix, the latter being singular
        ## (as it should be), as well as the equation of the hyperplane.

        dim(mcd$coeff) <- c(5, p)
        ans$cov <- ans$raw.cov <- mcd$initcovariance
        ans$center <- ans$raw.center <- as.vector(mcd$initmean)

        if(names && !is.null(nms <- dimn[[2]])) {
            dimnames(ans$cov) <- list(nms, nms)
            names(ans$center) <- nms
        }
        ans$n.obs <- n

	if(raw.only) {
	    ans$raw.only <- TRUE
	} else {
	    ## no longer relevant:
	    ##      if(mcd$exactfit == -1)
	    ##      stop("The program allows for at most ", mcd$kount, " observations.")
	    ##      if(mcd$exactfit == -2)
	    ##      stop("The program allows for at most ", mcd$kount, " variables.")
	    if(!(mcd$exactfit %in% c(1,2,3)))
		stop("Unexpected 'exactfit' code ", mcd$exactfit, ". Please report!")
	    ## new (2007-01) and *instead* of older long 'method' extension;
	    ## the old message is still *printed* via .MCDsingularityMsg()
	    ##
	    ## exactfit is now *passed* to result instead of coded into 'message':
	    ans$singularity <-
		list(kind = "on.hyperplane", exactCode = mcd$exactfit,
		     p = p, h = h, count = mcd$kount, coeff = mcd$coeff[1,])
	}
        ans$alpha <- alpha
        ans$quan <- h
        if(names && !is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$crit <- -Inf # = log(0)
        weights <- mcd$weights

      } ## end (raw.only || exact fit)

      else { ## have general position (exactfit == 0) : ------------------------

        ## FIXME? here, we assume that mcd$initcovariance is not singular:
        mah <- mahalanobis(x, mcd$initmean, mcd$initcovariance, tol = tolSolve)
        weights <- wgtFUN(mah)
        sum.w <- sum(weights)
        ans <- c(ans, cov.wt(x, wt = weights, cor=cor))
        ## simple check for singularity, much cheaper than determinant() below:
        sing.rewt <- any(apply(ans$cov == 0, 2, all))

        ## Compute and apply the consistency correction factor for
        ## the reweighted cov
        if(!sing.rewt && sum.w != n) {
	    cdelta.rew <- .MCDcons(p, sum.w/n) ## VT::19.3.2007
	    correct.rew <- if(use.correction) .MCDcnp2.rew(p, n, alpha) else 1.
	    cnp2 <- c(cdelta.rew, correct.rew)
	    ans$cov <- cdelta.rew * correct.rew * ans$cov
        }

        ##vt:: add also the best found subsample to the result list
        ans$best <- sort(as.vector(mcd$best))

        ans$alpha <- alpha
        ans$quan <- h
        ans$raw.cov <- mcd$initcovariance
        ans$raw.center <- as.vector(mcd$initmean)
        if(names && !is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$raw.weights <- weights
        ans$crit <- mcd$mcdestimate # now in log scale!
        ## 'mah' already computed above
        ans$raw.mah <- mah ## mahalanobis(x, ans$raw.center, ans$raw.cov, tol = tolSolve)
        ## Check if the reweighted scatter matrix is singular.
        if(sing.rewt || - determinant(ans$cov, logarithm = TRUE)$modulus[1]/p > logdet.Lrg) {
	    ans$singularity <- list(kind = paste0("reweighted.MCD",
				    if(sing.rewt)"(zero col.)"))
	    ans$mah <- mah
        }
        else {
            mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
            ans$mah <- mah
            weights <- wgtFUN(mah)
        }
      } ## end{ not exact fit }

    } ## end{ p >= 2 }

    ans$mcd.wt <- rep.int(NA, length(ok))
    ans$mcd.wt[ok] <- weights
    if(names) {
	if(length(dimn[[1]]))
	    names(ans$mcd.wt) <- dimn[[1]]
	if(length(dimn[[1]]))
	    dimnames(x)[[1]] <- names(ans$mcd.wt)[ok]
	else
	    dimnames(x) <- list(seq(along = ok)[ok], NULL)
    }
    ans$X <- x
    ans$wt <- NULL
    if(trace)
        cat(ans$method, "\n")
    ans$raw.cnp2 <- raw.cnp2
    ans$cnp2 <- cnp2
    if(nsamp == "deterministic")
	ans <- c(ans, mcd[c("iBest","n.csteps", if(save.hsets) "initHsets")])
    class(ans) <- "mcd"
    ## warn if we have a singularity:
    if(is.list(ans$singularity))
	warning(paste(strwrap(.MCDsingularityMsg(ans$singularity, ans$n.obs)), collapse="\n"),
		domain=NA)
    ## return
    ans
} ## {covMcd}

smoothWgt <- function(x, c, h) {
    ## currently drops all attributes, dim(), names(), etc
    ## maybe add 'keep.attributes = FALSE' (and pass to C)
    .Call(R_wgt_flex, x, c, h)
}

##' Martin Maechler's simple proposal for an *adaptive* cutoff
##' i.e., one which does *not* reject outliers in good samples asymptotically
.MCDadaptWgt.c <- function(n,p) {
    eps <- 0.4 / n ^ 0.6 # => 1-eps(n=100) ~= 0.975; 1-eps(n=10) ~= 0.90
    ## using upper tail:
    qchisq(eps, p, lower.tail=FALSE)
}


## Default wgtFUN()s :
.wgtFUN.covMcd <-
    list("01.original" = function(p, ...) {
	     cMah <- qchisq(0.975, p)
	     function(d) as.numeric(d < cMah)
	 },
	 "01.flex" = function(p, n, control) { ## 'control$beta' instead of 0.975
	     ## FIXME: update rrcov.control() to accept 'beta'
	     stopifnot(is.1num(beta <- control$beta), 0 <= beta, beta <= 1)
	     cMah <- qchisq(beta, p)
	     function(d) as.numeric(d < cMah)
	 },
	 "01.adaptive" = function(p, n, ...) { ## 'beta_n' instead of 0.975
	     cMah <- .MCDadaptWgt.c(n,p)
	     function(d) as.numeric(d < cMah)
	 },
	 "sm1.orig" = function(p, n, ...) {
	     cMah <- qchisq(0.975, p)
	     function(d) smoothWgt(d, c = cMah, h = 1)
	 },
	 "sm2.orig" = function(p, n, ...) {
	     cMah <- qchisq(0.975, p)
	     function(d) smoothWgt(d, c = cMah, h = 2)
	 },
	 "sm1.adaptive" = function(p, n, ...) {
	     cMah <- .MCDadaptWgt.c(n,p)
	     function(d) smoothWgt(d, c = cMah, h = 1)
	 },
	 "sm2.adaptive" = function(p, n, ...) {
	     cMah <- .MCDadaptWgt.c(n,p)
	     function(d) smoothWgt(d, c = cMah, h = 2)
	 },
	 "smE.adaptive" = function(p, n, ...) {
	     cMah <- .MCDadaptWgt.c(n,p)
	     ## TODO: find "theory" for h = f(cMah), or better c=f1(n,p); h=f2(n,p)
	     function(d) smoothWgt(d, c = cMah, h = max(2, cMah/4))
	 }
	 )


.MCDsingularityMsg <- function(singList, n.obs)
{
    stopifnot(is.list(singList))

    switch(singList$kind,
       "classical" = {
           "The classical covariance matrix is singular."
       },
       "reweighted.MCD" = {
           "The reweighted MCD scatter matrix is singular."
       },
       "identicalObs" = {
           sprintf("Initial scale 0 because more than 'h' (=%d) observations are identical.",
               singList$q)
       },
       "on.hyperplane" = {
           stopifnot(c("p", "count", "coeff") %in% names(singList))
	   obsMsg <- function(m, n)
	       paste("There are", m,
		     "observations (in the entire dataset of", n, "obs.)",
		     "lying on the")
           with(singList,
                c(switch(exactCode,
                         ## exactfit == 1 :
                         "The covariance matrix of the data is singular.",
                         ## exactfit == 2 :
                         c("The covariance matrix has become singular during",
                           "the iterations of the MCD algorithm."),
			 ## exactfit == 3:
			 paste0("The ", h,
				"-th order statistic of the absolute deviation of variable ",
				which(singList$coeff == 1), " is zero.")),

                  if(p == 2) {
                      paste(obsMsg(count, n.obs), "line with equation ",
                            signif(coeff[1], digits= 5), "(x_i1-m_1) +",
                            signif(coeff[2], digits= 5), "(x_i2-m_2) = 0",
                            "with (m_1,m_2) the mean of these observations.")
                  }
                  else if(p == 3) {
                      paste(obsMsg(count, n.obs), "plane with equation ",
                            signif(coeff[1], digits= 5), "(x_i1-m_1) +",
                            signif(coeff[2], digits= 5), "(x_i2-m_2) +",
                            signif(coeff[3], digits= 5), "(x_i3-m_3) = 0",
                            "with (m_1,m_2) the mean of these observations."
                            )
                  }
                  else { ##  p > 3 -----------
                      paste(obsMsg(count, n.obs), "hyperplane with equation ",
                                "a_1*(x_i1 - m_1) + ... + a_p*(x_ip - m_p) = 0",
                            "with (m_1, ..., m_p) the mean of these observations",
			    "and coefficients a_i from the vector   a <- ",
			    paste(deparse(zapsmall(coeff)), collapse="\n "))
                  }))
       },
       ## Otherwise
       stop("illegal 'singularity$kind'")
       ) ## end{switch}
}

nobs.mcd <- function (object, ...) object$n.obs

print.mcd <- function(x, digits = max(3, getOption("digits") - 3), print.gap = 2, ...)
{
    cat("Minimum Covariance Determinant (MCD) estimator approximation.\n",
        "Method: ", x$method, "\n", sep="")
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    if(is.list(x$singularity))
        cat(strwrap(.MCDsingularityMsg(x$singularity, x$n.obs)), sep ="\n")

    if(identical(x$nsamp, "deterministic"))
	cat("iBest: ", pasteK(x$iBest), "; C-step iterations: ", pasteK(x$n.csteps),
            "\n", sep="")
    ## VT::29.03.2007 - solve a conflict with fastmcd() in package robust -
    ##      also returning an object of class "mcd"
    xx <- NA
    if(!is.null(x$crit))
	xx <- format(x$crit, digits = digits)
    else if (!is.null(x$raw.objective))
	xx <- format(log(x$raw.objective), digits = digits)
    cat("Log(Det.): ", xx , "\n\nRobust Estimate of Location:\n")
    print(x$center, digits = digits, print.gap = print.gap, ...)
    cat("Robust Estimate of Covariance:\n")
    print(x$cov, digits = digits, print.gap = print.gap, ...)
    invisible(x)
}

summary.mcd <- function(object, ...)
{
    class(object) <- c("summary.mcd", class(object))
    object
}

print.summary.mcd <-
    function(x, digits = max(3, getOption("digits") - 3), print.gap = 2, ...)
{
    print.mcd(x, digits = digits, print.gap = print.gap, ...) # see above

    ## hmm, maybe not *such* a good idea :
    if(!is.null(x$cor)) {
        cat("\nRobust Estimate of Correlation: \n")
        dimnames(x$cor) <- dimnames(x$cov)
        print(x$cor, digits = digits, print.gap = print.gap, ...)
    }

    cat("\nEigenvalues:\n")
    print(eigen(x$cov, only.values = TRUE)$values, digits = digits, ...)

    if(!is.null(x$mah)) {
	cat("\nRobust Distances: \n")
	print(summary(x$mah, digits = digits), digits = digits, ...)
    }
    if(!is.null(wt <- x$mcd.wt))
	summarizeRobWeights(wt, digits = digits)
    invisible(x)
}

## NOTE:  plot.mcd() is in ./covPlot.R !
## ----                    ~~~~~~~~~~~

### Consistency and Finite Sample Correction Factors
###  .MCDcons()         .MCDcnp2() & .MCDcnp2.rew()

### now exported and documented in ../man/covMcd.Rd

##' Compute the consistency correction factor for the MCD estimate
##'    (see calfa in Croux and Haesbroeck)
##' @param p
##' @param alpha alpha ~= h/n = quan/n
##'    also use for the reweighted MCD, calling with alpha = 'sum(weights)/n'
MCDcons <- # <- *not* exported, but currently used in pkgs rrcov, rrcovNA
.MCDcons <- function(p, alpha)
{
    qalpha <- qchisq(alpha, p)
    caI <- pgamma(qalpha/2, p/2 + 1) / alpha
    1/caI
}

MCDcnp2 <- # <- *not* exported, but currently used in pkg rrcovNA
##' Finite sample correction factor for raw MCD:
.MCDcnp2 <- function(p, n, alpha)
{
    stopifnot(0 <= alpha, alpha <= 1, length(alpha) == 1)

    if(p > 2) {
	##				"alfaq"	        "betaq"	    "qwaarden"
        coeffqpkwad875 <- matrix(c(-0.455179464070565, 1.11192541278794, 2,
                                   -0.294241208320834, 1.09649329149811, 3),
                                 ncol = 2)
        coeffqpkwad500 <- matrix(c(-1.42764571687802,  1.26263336932151, 2,
                                   -1.06141115981725,  1.28907991440387, 3),
                                 ncol = 2)

        y.500 <- log( - coeffqpkwad500[1, ] / p^coeffqpkwad500[2, ] )
        y.875 <- log( - coeffqpkwad875[1, ] / p^coeffqpkwad875[2, ] )

        A.500 <- cbind(1, - log(coeffqpkwad500[3, ] * p^2))
        A.875 <- cbind(1, - log(coeffqpkwad875[3, ] * p^2))
        coeffic.500 <- solve(A.500, y.500)
        coeffic.875 <- solve(A.875, y.875)
        fp.500.n <- 1 - exp(coeffic.500[1]) / n^coeffic.500[2]
        fp.875.n <- 1 - exp(coeffic.875[1]) / n^coeffic.875[2]
    }
    else if(p == 2) {
        fp.500.n <- 1 - exp( 0.673292623522027) / n^0.691365864961895
        fp.875.n <- 1 - exp( 0.446537815635445) / n^1.06690782995919
    } else if(p == 1) {
        fp.500.n <- 1 - exp( 0.262024211897096) / n^0.604756680630497
        fp.875.n <- 1 - exp(-0.351584646688712) / n^1.01646567502486
    }

    ## VT:18.04.2007 - use simulated correction factors for several p and n:
    ## p in [1, 20] n in [2*p, ...]
    if(alpha == 0.5 && !is.na(fp.x <- MCDcnp2s$sim.0(p, n)))
        fp.500.n <- 1/fp.x

    fp.alpha.n <-
        if(alpha <= 0.875)
            fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
        else ##  0.875 < alpha <= 1
            fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)

    1/fp.alpha.n
} ## end{.MCDcnp2 }

MCDcnp2.rew <- # <- *not* exported, but currently used in pkg rrcovNA
##' Finite sample correction factor for *REW*eighted MCD
.MCDcnp2.rew <- function(p, n, alpha)
{
    stopifnot(0 <= alpha, alpha <= 1, length(alpha) == 1)

    if(p > 2) {
        ##                              "alfaq"         "betaq"        "qwaarden"
        coeffrewqpkwad875 <- matrix(c(-0.544482443573914, 1.25994483222292, 2,
                                      -0.343791072183285, 1.25159004257133, 3),
                                    ncol = 2)
        coeffrewqpkwad500 <- matrix(c(-1.02842572724793,  1.67659883081926, 2,
                                      -0.26800273450853,  1.35968562893582, 3),
                                    ncol = 2)

        y.500 <- log( - coeffrewqpkwad500[1, ] / p^ coeffrewqpkwad500[2, ] )
        y.875 <- log( - coeffrewqpkwad875[1, ] / p^ coeffrewqpkwad875[2, ] )

        A.500 <- cbind(1, - log(coeffrewqpkwad500[3, ] * p^2))
        coeffic.500 <- solve(A.500, y.500)
        A.875 <- cbind(1, - log(coeffrewqpkwad875[3, ] * p^2))
        coeffic.875 <- solve(A.875, y.875)
        fp.500.n <- 1 - exp(coeffic.500[1]) / n^ coeffic.500[2]
        fp.875.n <- 1 - exp(coeffic.875[1]) / n^ coeffic.875[2]
    }
    else if(p == 2) {
        fp.500.n <- 1 - exp( 3.11101712909049 ) / n^ 1.91401056721863
        fp.875.n <- 1 - exp( 0.79473550581058 ) / n^ 1.10081930350091
    } else if(p == 1) {
        fp.500.n <- 1 - exp( 1.11098143415027 ) / n^ 1.5182890270453
        fp.875.n <- 1 - exp( -0.66046776772861) / n^ 0.88939595831888
    }

    ## VT:18.04.2007 - use simulated correction factors for several p and n:
    ## p in [1, 20] n in [2*p, ...]
    if(alpha == 0.5 && !is.na(fp.x <- MCDcnp2s$sim.rew(p, n)))
        fp.500.n <- 1/fp.x

    fp.alpha.n <-
        if(alpha <= 0.875)
            fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
        else ##  0.875 < alpha <= 1
            fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)

    1/fp.alpha.n
} ## end{.MCDcnp2.rew }


.fastmcd <- function(x, h, nsamp, nmini, kmini, trace = 0)
{
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

    ##   parameters for partitioning {equal to those in Fortran !!}
    ## kmini <- 5
    ## nmini <- 300
    stopifnot(length(kmini <- as.integer(kmini)) == 1, kmini >= 2L,
              is.1num(nmini), is.finite(nmaxi <- as.double(nmini)*kmini),
              nmaxi * p < .Machine$integer.max)
    nmaxi <- as.integer(nmaxi)
    km10 <- 10*kmini

    ## vt::03.02.2006 - added options "best" and "exact" for nsamp
    ##
    nLarge <- 100000 # was 5000 before Nov.2009 -- keep forever now; user can say "exact"
    if(is.numeric(nsamp) && (nsamp < 0 || nsamp == 0 && p > 1)) {
        nsamp <- -1
    } else if(nsamp == "exact" || nsamp == "best") {
        if(n > 2*nmini-1) {
            warning("Options 'best' and 'exact' not allowed for n greater than  2*nmini-1 =",
                    2*nmini-1,".\nUsing default.\n")
            nsamp <- -1
        } else {
	    myk <- p + 1 ## was 'p'; but p+1 ("nsel = nvar+1") is correct
            nall <- choose(n, myk)
            msg <- paste("subsets of size", myk, "out of", n)
            if(nall > nLarge && nsamp == "best") {
                nsamp <- nLarge
                warning("'nsamp = \"best\"' allows maximally ",
                        format(nLarge, scientific=FALSE),
                        " subsets;\ncomputing these ", msg,
                        immediate. = TRUE)
            } else {   ## "exact" or ("best"  &  nall < nLarge)
                nsamp <- 0 ## all subsamples -> special treatment in Fortran
                if(nall > nLarge) {
                    msg <- paste("Computing all", nall, msg)
                    if(nall > 10*nLarge)
                        warning(msg, "\n This may take a",
                                if(nall/nLarge > 100) " very", " long time!\n",
                                immediate. = TRUE)
                    else message(msg)
                }
            }
        }
    }

    if(!is.numeric(nsamp) || nsamp == -1) { # still not defined
        ## set it to the default :
        nsamp.def <- rrcov.control()$nsamp
	warning(gettextf("Invalid number of trials nsamp=%s. Using default nsamp=%d.",
			 format(nsamp), nsamp.def),
		domain=NA)
        nsamp <- nsamp.def
    }

    if(nsamp > (mx <- .Machine$integer.max)) {
	warning("nsamp > i_max := maximal integer -- not allowed;\n",
		" set to i_max = ", mx)
	nsamp <- mx
    }

    ##   Allocate temporary storage for the Fortran implementation,
    ##   directly in the .Fortran() call.
    ##    (if we used C + .Call() we would allocate all there, and be quite faster!)

    .Fortran(rffastmcd,
             x = if(is.double(x)) x else as.double(x),
             n =    as.integer(n),
             p =    as.integer(p), ## = 'nvar'  in Fortran
             nhalff =   as.integer(h),
             nsamp  =   as.integer(nsamp), # = 'krep'
             nmini  =   as.integer(nmini),
	     kmini  =	kmini,
             initcovariance = double(p * p),
             initmean       = double(p),
             best       = rep.int(as.integer(10000), h),
             mcdestimate = double(1), ## = 'det'
             weights   = integer(n),
             exactfit  = integer(1), # output indicator: 0: ok; 1: ..., 2: ..
             coeff     = matrix(double(5 * p), nrow = 5, ncol = p), ## plane
             kount     = integer(1),
             adjustcov = double(p * p), ## used in ltsReg() !
             ## integer(1), ## << 'seed' no longer used
             temp   = integer(n),
             index1 = integer(n),
             index2 = integer(n),
             indexx = integer(n),
             nmahad = double(n),
             ndist  = double(n),
             am     = double(n),
             am2    = double(n),
             slutn  = double(n),

             med   = double(p),
             mad   = double(p),
             sd    = double(p),
             means = double(p),
             bmeans= double(p),
             w     = double(p),
             fv1   = double(p),
             fv2   = double(p),

             rec   = double(p+1),
             sscp1 = double((p+1)*(p+1)),
             cova1 = double(p * p),
             corr1 = double(p * p),
             cinv1 = double(p * p),
             cova2 = double(p * p),
             cinv2 = double(p * p),
             z     = double(p * p),

             cstock = double(10 * p * p), # (10,nvmax2)
             mstock = double(10 * p),     # (10,nvmax)
             c1stock = double(km10 * p * p), # (km10,nvmax2)
             m1stock = double(km10 * p),     # (km10,nvmax)
             dath = double(nmaxi * p),       # (nmaxi,nvmax)

             cutoff = qchisq(0.975, p),
             chimed = qchisq(0.5,   p),
             i.trace= as.integer(trace)
             )[ ## keep the following ones:
               c("initcovariance", "initmean", "best", "mcdestimate",
                 "weights", "exactfit", "coeff", "kount", "adjustcov") ]
}

##
## VT:18.04.2007 - use simulated correction factors for several p and n
##   and alpha = 1/2  (the default in rrcov.control())
##       ~~~~~~~~~~~
##  p in [1, 20] n in [2*p, ...]
##  see the modifications in.MCDcnp2() and.MCDcnp2.rew
##

##  VT::08.06.2007 - fixed the simulated values (especially for p=1)
##  VT::11.05.2007 - reduce the usage of the simulated correction factors to only those that
##  are definitvely wrong (negative or very large). This is done by:
##      a) reducing p.max
##      b) reducing n.max
##  NB: In general, "wrong" are the factors for the reweighted matrix, but whenever a simulated
##      value for the reweighted is used, the corresponding simulated must be used for the raw too.
##

## MM::2014-04 :
MCDcnp2s <- local({
    p.min <- 1L
    p.max <- 9L # was 20
    ncol <- 20L # the number of column in the matrices
    n.min <- as.integer(
###  p =  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
        c(1,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40))
    n.max <- as.integer(
        c(2,  6, 10, 13, 16, 18, 20, 20, 20, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60))
##was  c(22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60)
## these are the right (simulated) values for n.max

    n.min.rew <- n.min
    n.max.rew <- n.max

    m.0 <- matrix(
        c(1, 3.075819, 1.515999, 2.156169, 1.480742, 1.765485, 1.460206, 1.603707, 1.427429, 1.504712, 1.334528, 1.48297,  1.355308, 1.383867, 1.319241, 1.36065,  1.307467, 1.365596, 1.255259, 1.352741, 1.239381, 3.15342, 1.799889, 2.258497, 1.688312, 1.906779, 1.548203, 1.724785, 1.500873, 1.573442, 1.417137, 1.540805, 1.395945, 1.472596, 1.394247, 1.377487, 1.337394, 1.369354, 1.333378, 1.3181, 1.313813, 1.315528, 2.12777, 2.718898, 1.993509, 2.220433, 1.820585, 1.97782, 1.672455, 1.770151, 1.587478, 1.685352, 1.539295, 1.584536, 1.499487, 1.50702, 1.41952, 1.449058, 1.393042, 1.432999, 1.369964, 1.400997, 1.333824, 2.950549, 2.145387, 2.382224, 1.927077, 2.032489, 1.8371, 1.877833, 1.710891, 1.756053, 1.620778, 1.657761, 1.558978, 1.56257, 1.508633, 1.534406, 1.46709, 1.468734, 1.432529, 1.455283, 1.386975, 1.417532, 2.229573, 2.494447, 2.016117, 2.190061, 1.877996, 1.978964, 1.767284, 1.836948, 1.677372, 1.743316, 1.616383, 1.655964, 1.55484, 1.594831, 1.502185, 1.543723, 1.467005, 1.491123, 1.44402, 1.446915, 1.401578, 2.580264, 2.109121, 2.240741, 1.944719, 2.043397, 1.821808, 1.89725, 1.748788, 1.786988, 1.659333, 1.697012, 1.610622, 1.616503, 1.538529, 1.562024, 1.499964, 1.529344, 1.474519, 1.483264, 1.441552, 1.434448, 2.165233, 2.320281, 2.007836, 2.086471, 1.884052, 1.950563, 1.76926, 1.843328, 1.708941, 1.741039, 1.627206, 1.644755, 1.580563, 1.593402, 1.527312, 1.568418, 1.501462, 1.502542, 1.464583, 1.467921, 1.431141, 2.340443, 2.048262, 2.161097, 1.926082, 1.995422, 1.81446, 1.853165, 1.738533, 1.784456, 1.679444, 1.696463, 1.612931, 1.629483, 1.548186, 1.580026, 1.52198, 1.531111, 1.482914, 1.484824, 1.442726, 1.447838, 2.093386, 2.185793, 1.948989, 2.02804, 1.867137, 1.907732, 1.771923, 1.800413, 1.691612, 1.720603, 1.642705, 1.649769, 1.589028, 1.598955, 1.539759, 1.55096, 1.503965, 1.50703, 1.471349, 1.469791, 1.436959, 2.218315, 1.997369, 2.041128, 1.887059, 1.928524, 1.79626, 1.827538, 1.716748, 1.735696, 1.658329, 1.664211, 1.599286, 1.611511, 1.553925, 1.562637, 1.516805, 1.529894, 1.476064, 1.482474, 1.453253, 1.458467, 2.0247, 2.07899, 1.921976, 1.949376, 1.824629, 1.851671, 1.744713, 1.765647, 1.683525, 1.685592, 1.625113, 1.624961, 1.571921, 1.581223, 1.535257, 1.537464, 1.497165, 1.504879, 1.468682, 1.469319, 1.448344, 2.092315, 1.941412, 1.969843, 1.844093, 1.866133, 1.766145, 1.783829, 1.703613, 1.709714, 1.646078, 1.654264, 1.594523, 1.598488, 1.545105, 1.555356, 1.514627, 1.521353, 1.483958, 1.487677, 1.449191, 1.459721, 1.958987, 1.985144, 1.87739, 1.879643, 1.786823, 1.799642, 1.720015, 1.724688, 1.663539, 1.662997, 1.609267, 1.615124, 1.56746, 1.562026, 1.520586, 1.52503, 1.493008, 1.502496, 1.471983, 1.468546, 1.435064, 1.994706, 1.880348, 1.894254, 1.805827, 1.815965, 1.744296, 1.743389, 1.665481, 1.681644, 1.624466, 1.626109, 1.584028, 1.5818, 1.54376, 1.547237, 1.504878, 1.515087, 1.479032, 1.47936, 1.450758, 1.45073, 1.892685, 1.91087, 1.825301, 1.827176, 1.745363, 1.746115, 1.693373, 1.701692, 1.648247, 1.637112, 1.594648, 1.592013, 1.554849, 1.55013, 1.522186, 1.520901, 1.492606, 1.493072, 1.460868, 1.46733, 1.440956, 1.92771, 1.835696, 1.841979, 1.775991, 1.766092, 1.703807, 1.708791, 1.654985, 1.655917, 1.602388, 1.611867, 1.570765, 1.573368, 1.53419, 1.529033, 1.506767, 1.503596, 1.481126, 1.471806, 1.444917, 1.451682, 1.850262, 1.855034, 1.778997, 1.789995, 1.718871, 1.717326, 1.667357, 1.666291, 1.619743, 1.631475, 1.582624, 1.58766, 1.546302, 1.545063, 1.512222, 1.517888, 1.489127, 1.487271, 1.466722, 1.463618, 1.444137, 1.8709, 1.794033, 1.80121, 1.736376, 1.740201, 1.673776, 1.682541, 1.638153, 1.642294, 1.604417, 1.597721, 1.559534, 1.559108, 1.533942, 1.529348, 1.499517, 1.501586, 1.473147, 1.473031, 1.457615, 1.452348, 1.805753, 1.812952, 1.746549, 1.747222, 1.696924, 1.694957, 1.652157, 1.650568, 1.607807, 1.613666, 1.577295, 1.570712, 1.543704, 1.538272, 1.515369, 1.517113, 1.487451, 1.491593, 1.464514, 1.464658, 1.439359, 1.823222, 1.758781, 1.767358, 1.70872, 1.712926, 1.666956, 1.667838, 1.62077, 1.621445, 1.592891, 1.58549, 1.55603, 1.559042, 1.521501, 1.523342, 1.499913, 1.501937, 1.473359, 1.472522, 1.452613, 1.452448),
                         ncol = ncol)

    m.rew <- matrix(
    c(1, 0.984724, 0.970109, 0.978037, 0.979202, 0.982933, 1.001461, 1.026651, 0.981233, 1.011895, 1.017499, 0.964323, 1.026574, 1.006594, 0.980194, 1.009828, 0.998083, 0.966173, 1.009942, 0.99916, 1.021521, 2.216302, 1.418526, 1.635601, 1.31402, 1.33975, 1.251798, 1.210917, 1.133114, 1.150666, 1.138732, 1.096822, 1.076489, 1.058343, 1.045746, 1.036743, 1.008929, 1.049537, 1.028148, 1.027297, 1.020578, 1.00074, 1.73511, 2.06681, 1.545905, 1.659655, 1.456835, 1.47809, 1.331966, 1.334229, 1.231218, 1.220443, 1.198143, 1.193965, 1.142156, 1.146231, 1.124661, 1.112719, 1.089973, 1.070606, 1.082681, 1.061243, 1.053191, 2.388892, 1.847626, 1.96998, 1.630723, 1.701272, 1.521008, 1.553057, 1.382168, 1.414555, 1.326982, 1.321403, 1.265207, 1.264856, 1.200418, 1.21152, 1.17531, 1.168536, 1.140586, 1.14457, 1.111392, 1.112031, 1.968153, 2.168931, 1.784373, 1.894409, 1.667912, 1.693007, 1.545176, 1.582428, 1.45319, 1.480559, 1.371611, 1.358541, 1.330235, 1.30264, 1.257518, 1.244156, 1.221907, 1.22455, 1.178965, 1.177855, 1.166319, 2.275891, 1.866587, 2.014249, 1.750567, 1.829363, 1.650019, 1.689043, 1.562539, 1.561359, 1.473378, 1.488554, 1.411097, 1.416527, 1.35117, 1.361044, 1.30205, 1.299037, 1.250265, 1.260083, 1.218665, 1.236027, 1.95771, 2.074066, 1.847385, 1.905408, 1.71393, 1.768425, 1.63908, 1.67234, 1.564992, 1.562337, 1.49229, 1.499573, 1.420813, 1.424067, 1.383947, 1.378726, 1.33062, 1.330071, 1.279404, 1.295302, 1.263947, 2.164121, 1.871024, 1.979485, 1.782417, 1.84489, 1.706023, 1.734857, 1.622782, 1.634869, 1.55196, 1.554423, 1.482325, 1.509195, 1.440726, 1.436328, 1.386335, 1.396277, 1.347939, 1.346732, 1.310242, 1.309371, 1.938822, 2.050409, 1.834863, 1.882536, 1.737494, 1.761608, 1.65742, 1.687579, 1.591863, 1.60158, 1.520982, 1.535234, 1.470649, 1.486485, 1.42892, 1.435574, 1.384132, 1.382329, 1.343281, 1.346581, 1.315111, 2.063894, 1.880094, 1.907246, 1.78278, 1.806648, 1.6952, 1.720922, 1.63084, 1.635274, 1.565423, 1.56171, 1.512015, 1.4986, 1.463903, 1.456588, 1.422856, 1.407325, 1.376724, 1.373923, 1.346464, 1.34259, 1.898389, 1.950406, 1.812053, 1.849175, 1.72649, 1.737651, 1.646719, 1.655112, 1.587601, 1.597894, 1.539877, 1.53329, 1.495054, 1.490548, 1.445249, 1.446037, 1.410272, 1.412274, 1.375797, 1.369604, 1.341232, 1.992488, 1.830452, 1.857314, 1.758686, 1.763822, 1.683215, 1.679543, 1.619269, 1.608512, 1.565, 1.562282, 1.498869, 1.51325, 1.470912, 1.464654, 1.427573, 1.439301, 1.402308, 1.391006, 1.37074, 1.367573, 1.855502, 1.891242, 1.77513, 1.790618, 1.706443, 1.713098, 1.642896, 1.636577, 1.580366, 1.581752, 1.542937, 1.531668, 1.487894, 1.492039, 1.460304, 1.449762, 1.4219, 1.420953, 1.390137, 1.388677, 1.360506, 1.908277, 1.802091, 1.806128, 1.723757, 1.727249, 1.659883, 1.670056, 1.605209, 1.611481, 1.558846, 1.551762, 1.512951, 1.511515, 1.468948, 1.476073, 1.441508, 1.434997, 1.412687, 1.406782, 1.380452, 1.375924, 1.811415, 1.822311, 1.740544, 1.739355, 1.68127, 1.685342, 1.620281, 1.622572, 1.579611, 1.570103, 1.529881, 1.530097, 1.490041, 1.4947, 1.457329, 1.456344, 1.423363, 1.428653, 1.399988, 1.390069, 1.376594, 1.837723, 1.76039, 1.771031, 1.697404, 1.690915, 1.634409, 1.63713, 1.589594, 1.586521, 1.552974, 1.545571, 1.505923, 1.512794, 1.477833, 1.477821, 1.444241, 1.44452, 1.419258, 1.421297, 1.394924, 1.389393, 1.779716, 1.781271, 1.706031, 1.71224, 1.655099, 1.654284, 1.608878, 1.605955, 1.565683, 1.565938, 1.523594, 1.531235, 1.492749, 1.486786, 1.457635, 1.461416, 1.432472, 1.430164, 1.404441, 1.400021, 1.378273, 1.798932, 1.735577, 1.727031, 1.671049, 1.677601, 1.624427, 1.617626, 1.579533, 1.579987, 1.544635, 1.538715, 1.504538, 1.50726, 1.477163, 1.477084, 1.450861, 1.444496, 1.428416, 1.422813, 1.400185, 1.39552, 1.750193, 1.752145, 1.690365, 1.692051, 1.642391, 1.63858, 1.600144, 1.596401, 1.558305, 1.555932, 1.525968, 1.522984, 1.491563, 1.492554, 1.467575, 1.45786, 1.437545, 1.430893, 1.413983, 1.409386, 1.391943, 1.762922, 1.701346, 1.704996, 1.6556, 1.655548, 1.611964, 1.615219, 1.569103, 1.571079, 1.540617, 1.541602, 1.503791, 1.50195, 1.478069, 1.47678, 1.452458, 1.451732, 1.429144, 1.426547, 1.40363, 1.402647),
                         ncol = ncol)

    rm(ncol)
    list(
        sim.0 = function(p, n) {
            p. <- p - p.min + 1L
            if(p.min     <= p && p <= p.max &&
               n.min[p.] <= n && n <= n.max[p.]) {
                nind <- n - n.min[p.] + 1L
                m.0[nind, p.]
                ##=
            } else NA
        },
        sim.rew = function(p, n) {
            p. <- p - p.min + 1L
            if(p.min         <= p && p <= p.max &&
               n.min.rew[p.] <= n && n <= n.max.rew[p.]) {
                nind <- n - n.min.rew[p.] + 1L
                m.rew[nind, p.]
                ##===
            } else NA
        })
}) ## end{MCDcnp2s}

if(FALSE) { ## For experimentation:
    ls.str( ee <- environment(MCDcnp2s$sim.0) )
    matplot(1:21, ee$m.0, type = "o", xlab = "p - p.min + 1")
}
