### -*- mode: R ; delete-old-versions: never -*-

##' Computes the MCD estimator of a multivariate data sets in a \emph{deterministic}
##' way.
##'
##' The MCD estimator is given by the subset of h observations with smallest
##' covariance determinant.  The MCD location estimate is then
##' the mean of those h points, and the MCD scatter estimate is
##' their covariance matrix.  The default value of h is roughly
##' 0.75n (where n is the total number of observations), but the
##' user may choose each value between n/2 and n. Based on the
##' raw estimates, weights are assigned to the observations such
##' that outliers get zero weight. The reweighted MCD estimator
##' is then given by the mean and covariance matrix of the cases
##' with non-zero weight.
##
##' To compute an approximate MCD estimator deterministically, six initial robust h-subsets are
##'  constructed based on robust transformations of variables or robust and
##'  fast-to-compute estimators of multivariate location and shape. Then
##'  C-steps are applied on these h-subsets until convergence. Note that the
##'  resulting algorithm is not fully affine equivariant, but it is often
##'  faster than the FAST-MCD algorithm which is affine equivariant
##'  (see covMcd()).
##'  Note that this function can not handle exact fit situations: if the
##'  raw covariance matrix is singular, the program is stopped. In that
##'  case, it is recommended to apply the covMcd() function.
##'
##'  The MCD method is intended for continuous variables, and assumes that
##'  the number of observations n is at least 5 times the number of variables p.
##'  If p is too large relative to n, it would be better to first reduce
##'  p by variable selection or robust principal components (see the functions
##'  robust principal components in package 'rrcov').
##'
##' @title Compute the MCD estimator of multivariate data in a deterministic way
##' @references
##'  Hubert, M., Rousseeuw, P.J. and Verdonck, T. (2012),
##'  "A deterministic algorithm for robust location and scatter", Journal of
##'  Computational and Graphical Statistics, in press.
##' @param x a numerical matrix. The columns represent variables, and rows represent observations.
##' @param h The quantile of observations whose covariance determinant will
##'          be minimized.  Any value between n/2 and n may be specified.
##' @param hsets.init If one gives here already a matrix with for each column an
##' ordering of the observations (first the one with smallest statistical
##' distance), then the initial shape estimates are not calculated.
##' Default value = NULL.
##' @param save.hsets
##' @param full.h
##' @param scalefn function (or "rule") to estimate the scale.
##' @param maxcsteps
##' @param warn.nonconv.csteps
##' @param warn.wrong.obj.conv
##' @param trace
##' @param names
##' @return
##' @author Valentin Todorov;  many tweaks by Martin Maechler
.detmcd <- function(x, h, hsets.init=NULL,
		    save.hsets = missing(hsets.init), full.h = save.hsets,
		    scalefn, maxcsteps = 200,
		    warn.nonconv.csteps = getOption("robustbase:warn.nonconv.csteps", TRUE),
		    warn.wrong.obj.conv = getOption("robustbase:warn.wrong.obj.conv",FALSE),
		    trace = as.integer(trace), names = TRUE)
{
    stopifnot(length(dx <- dim(x)) == 2, h == as.integer(h), h >= 1)
    n <- dx[1]
    p <- dx[2]
    stopifnot(p >= 1, n >= 1)
    scalefn <- robScalefn(scalefn, n)

    ## kmini <- 5      # number of sub-data sets(if we use them some day)
    ##                 # for now we use it as number of rows in the returned
    ##                 # matrix 'coeff' for exact fit (also not used currently).

    ## cutoff <- qchisq(0.975, p)
    ## chimed <- qchisq(0.5,   p)

    ## Center and scale the data
    vnms <- colnames(x) # speedup only: store and put back at end
    z <- doScale(unname(x), center=median, scale=scalefn)
    z.center <- z$center
    z.scale <- z$scale
    z <- z$x

    ## Assume that 'hsets.init' already contains h-subsets: the first h observations each
    if(is.null(hsets.init)) {
	hsets.init <- r6pack(z, h=h, full.h=full.h, scaled=TRUE, scalefn=scalefn)
	dh <- dim(hsets.init)
    } else { ## user specified, (even just *one* vector):
	if(is.vector(hsets.init)) hsets.init <- as.matrix(hsets.init)
	dh <- dim(hsets.init)
	if(dh[1] < h || dh[2] < 1)
	    stop("'hsets.init' must be a  h' x L  matrix (h' >= h) of observation indices")
        ## TODO?: We could *extend* the sets to large h, even all n
        ## ====> could input the 'best' sets, also e.g. from fastmcd
	if(full.h && dh[1] != n)
	    warning("'full.h' is true, but 'hsets.init' has less than n rows")
##	    stop("When 'full.h' is true, user specified 'hsets.init' must have n rows")
	if(min(hsets.init) < 1 || max(hsets.init) > n)
	   stop("'hsets.init' must be in {1,2,...,n}; n = ", n)
    }

    nsets <- ncol(hsets.init)# typically 6, currently
    ## Some initializations.
    hset.csteps <- integer(nsets)
    bestobj <- Inf
    for(i in 1:nsets)
    {
        if(trace) {
            if(trace >= 2)
                cat(sprintf("H-subset %d = observations c(%s):\n-----------\n",
			    i, pasteK(hsets.init[1:h,i])))
            else
                cat(sprintf("H-subset %d: ", i))
        }
        for(j in 1:maxcsteps)
        {
            if(j == 1)  {
                obs_in_set <- hsets.init[1:h,i]     # start with the i-th initial set
            } else { # now using 'svd' from last step
		score <- (z - rep(svd$center, each=n)) %*% svd$loadings
                mah <- mahalanobisD(score, center=FALSE, sd = sqrt(abs(svd$eigenvalues)))
                obs_in_set <- sort.list(mah)[1:h] #, partial = 1:h not yet
            }
            ## [P,T,L,r,centerX,meanvct] = classSVD(data(obs_in_set,:));
	    svd <- classPC(z[obs_in_set, ,drop=FALSE], signflip=FALSE)
            obj <- sum(log(svd$eigenvalues))

	    if(svd$rank < p) { ## FIXME --> return exact fit property rather than stop() ??
		stop('More than h of the observations lie on a hyperplane.')
                ## TODO
                exactfit <- TRUE
                ## coeff <- ...
            }
            if(j >= 2 && obj == prevdet) {
                ## MM:: 2014-10-25: objective function check is *not* good enough:
                if(identical(obs_in_set, prevobs))
                    break
		## else :
		if(warn.wrong.obj.conv)
		    warning(sprintf(
			"original detmcd() wrongly declared c-step convergence (obj=%g, i=%d, j=%d)",
			obj, i,j))
            }
            prevdet <- obj
            prevobs <- obs_in_set
        }
        hset.csteps[i] <- j # how many csteps necessary to converge.
        if(trace) cat(sprintf("%3d csteps, obj=log(det|.|)=%g", j, obj))

        if(obj < bestobj) {
            if(trace) cat(" = new optim.\n")
            ## bestset           : the best subset for the whole data.
            ## bestobj           : objective value for this set.
            ## initmean, initcov : the mean and covariance matrix of this set
            bestset <- obs_in_set
            bestobj <- obj
            initmean <- svd$center
            L <- svd$loadings
            ## MM speedup: was L Diag L' = L %*% diag(svd$eigenvalues) %*% t(L)
	    initcov <- tcrossprod(L * rep(svd$eigenvalues, each=nrow(L)), L)
            ## raw.initcov <- initcov
            ## rew.Hsubsets.Hopt <- bestset
            ind.best <- i # to determine which subset gives best results.
	} else if(obj == bestobj) ## store as well:
	    ind.best <- c(ind.best, i)
        else if(trace) cat("\n")

    } ## for(i in 1:nsets)

    if(warn.nonconv.csteps && any(eq <- hset.csteps == maxcsteps)) {
	p1 <- paste(ngettext(sum(eq), "Initial set", "Initial sets"),
		    pasteK(which(eq)))
	warning(sprintf("%s did not converge in maxcsteps=%d concentration steps",
			p1, maxcsteps), domain=NA)
    }
    ## reweighting <- FALSE # it happens in covMcd()
    ## if(reweighting) {
    ##     svd <- classPC(z[bestset, ], signflip=FALSE)    # [P,T,L,r,centerX,meanvct] = classSVD(data(bestset,:));
    ##     mah <- mahalanobisD((z - rep(svd$center, each=n)) %*% svd$loadings,
    ##                        FALSE, sqrt(abs(svd$eigenvalues)))
    ##     sortmah <- sort(mah)
    ## }
    ##    factor <- sortmah[h]/qchisq(h/n, p)
    ##    raw.cov <- factor*initcov
    ## raw.cov <- initcov

    ## We express the results in the original units [restoring var.names]:
    raw.cov <- initcov * tcrossprod(z.scale)
    raw.center <- initmean * z.scale + z.center
    if(names) {
	dimnames(raw.cov) <- list(vnms, vnms)
	names(raw.center) <- vnms
    }
    raw.objective <- bestobj + 2*sum(log(z.scale)) # log(det = obj.best * prod(z.scale)^2)
    ## raw.mah <- mahalanobis(x, raw.center, raw.cov, tol=1E-14)
    ## medi2 <- median(raw.mah)

    list(initcovariance = raw.cov,
	 initmean = raw.center,
	 best = bestset,
	 mcdestimate = raw.objective, # determinant (goes to crit)
	 ## , weights=NULL# FIXME - goes to raw.weights
         iBest = ind.best,
         n.csteps = hset.csteps,
         initHsets = if(save.hsets) hsets.init,
         exactfit =  0 # <- FIXME
	 ## once we'd test for exact fit, we'd return:
	 ## , coeff=matrix(rep(0, kmini*p), nrow=kmini)
	 ## , kount=0 # FIXME
	 )
} ## .detmcd()

robScalefn <- function(scalefn, n) {
    if(missing(scalefn) || is.null(scalefn))
	scalefn <- .scalefn.default
    if(is.function(scalefn))
	scalefn
    else
	switch(scalefn,
	       ## Hubert, Rousseeuw, Verdonck, JCGS 2012 :
	       "hrv2012" = if(n < 1000) Qn else scaleTau2,
	       ## Version of 2014:
	       "v2014"	 = if(n < 5000) Qn else scaleTau2,
	       ## otherwise
	       stop(gettextf("Invalid  scalefn='%s': must be function or a valid string",
			     scalefn),
		    domain=NA))
}

doScale <- function (x, center, scale)
{
    stopifnot(is.numeric(p <- ncol(x)))
    ## MM: follow standard R's	scale.default() as much as possible

    centerFn <- is.function(center)
    doIt <- if(centerFn) {
        centerName <- deparse(substitute(center)) # "median" typically
	center <- apply(x, 2L, center)
	TRUE
    } else {
	if(length(center) == p && is.numeric(center))
	    TRUE
	else if(missing(center) || is.null(center)) {
            center <- 0; FALSE
	} else
	    stop(gettextf("'%s' must be a function, numeric vector of length p, or NULL",
                          "center"), domain=NA)
    }
    if(doIt)
	x <- sweep(x, 2L, center, `-`, check.margin=FALSE)

    scaleFn <- is.function(scale)
    doIt <- if(scaleFn) {
	scale <- apply(x, 2L, scale)
	TRUE
    } else {
	if(length(scale) == p && is.numeric(scale))
	    TRUE
	else if(missing(scale) || is.null(scale)) {
	    scale <- 1
	    FALSE
	} else
	    stop(gettextf("'%s' must be a function, numeric vector of length p, or NULL",
                          "scale"), domain=NA)
    }
    if(doIt) {
        if(any(is.na(scale)) || any(scale < 0))
            stop("provide better scale; must be all positive")
        if(any(s0 <- scale == 0)) {
## FIXME:
### Better and easier alternative (and as "FAST MCD"): return "singular cov.matrix"
### since scale 0 ==>  more than 50% points are on hyperplane x[,j] == const.
            ## find scale if there is any variation; otherwise use s := 1
            S <- if(centerFn && centerName == "median")
                     abs else function(.) abs(. - median(.))
            non0Q <- function(u) {
                alph <- c(10:19, 19.75)/20 # not all the way to '1' {=> finite qnorm()}
                qq <- quantile(S(u), probs=alph, names=FALSE)
                if(any(pos <- qq != 0)) { ## the first non-0 if there is one
                    i <- which.max(pos)
                    qq[i] / qnorm((alph[i] + 1)/2)
                } else 1
            }
            scale[s0] <- apply(x[,s0, drop=FALSE], 2L, non0Q)
        }
        x <- sweep(x, 2L, scale, `/`, check.margin = FALSE)
    }
    ## return
    list(x=x, center=center, scale=scale)
}

##' @title Robust Distance based observation orderings based on robust "Six pack"
##' @param x  n x p data matrix
##' @param h  integer
##' @param full.h full (length n) ordering or only the first h?
##' @param scaled is 'x' is already scaled?  otherwise, apply doScale(x, median, scalefn)
##' @param scalefn function to compute a robust univariate scale.
##' @return a h' x 6 matrix of indices from 1:n; if(full.h) h' = n else h' = h
r6pack <- function(x, h, full.h, scaled=TRUE,
                   scalefn = rrcov.control()$scalefn)
{
    ## As the considered initial estimators Sk may have very
    ## inaccurate eigenvalues, we try to 'improve' them by applying
    ## a transformation similar to that used in the OGK algorithm.
    ##
    ## After that compute the corresponding distances, order them and
    ## return the indices
    initset <- function(data, scalefn, P, h)
    {
        stopifnot(length(d <- dim(data)) == 2, length(h) == 1, h >= 1)
        n <- d[1]
        stopifnot(h <= n)
        lambda <- doScale(data %*% P, center=median, scale=scalefn)$scale
        sqrtcov    <- P %*% (lambda * t(P)) ## == P %*% diag(lambda) %*% t(P)
        sqrtinvcov <- P %*% (t(P) / lambda) ## == P %*% diag(1/lambda) %*% t(P)
	estloc <- colMedians(data %*% sqrtinvcov) %*% sqrtcov
        centeredx <- (data - rep(estloc, each=n)) %*% P
	sort.list(mahalanobisD(centeredx, FALSE, lambda))[1:h]# , partial = 1:h
    }

    ##
    ##  Compute the raw OGK estimator. For m(.) and s(.) (robust
    ##  univariate estimators of location and scale) use the median
    ##  and Qn for reasons of simplicity (no choice of tuning parameters)
    ##  and to be consistent with the other components of DetMCD.
    ##
    ogkscatter <- function(Y, scalefn, only.P = TRUE)
    {
        stopifnot(length(p <- ncol(Y)) == 1, p >= 1)
        U <- diag(p)

        for(i in seq_len(p)[-1L]) {# i = 2:p
            sYi <- Y[,i]
            ii <- seq_len(i - 1L)
            for(j in ii) {
                sYj <- Y[,j]
                U[i,j] <- (scalefn(sYi + sYj)^2 - scalefn(sYi - sYj)^2) / 4
            }
            ## also set the upper triangle
            U[ii,i] <- U[i,ii]
        }

        ## now done above: U <- lower.tri(U) * U + t(U)    #    U <- tril(U, -1) + t(U)
        P <- eigen(U, symmetric=TRUE)$vectors
	if(only.P)
	    return(P)
        ## else :
        Z <- Y %*% t(P)
        sigz <- apply(Z, 2, scalefn)
        lambda <- diag(sigz^2)

        list(P=P, lambda=lambda)
    }

    stopifnot(length(dx <- dim(x)) == 2)
    n <- dx[1]
    p <- dx[2]

    ## If scalefn is missing or is NULL, use Qn for smaller data sets (n < 1000)
    ## and tau-scale of Yohai and Zamar (1988) otherwise.
    scalefn <- robScalefn(scalefn, n)

    ## If the data was not scaled already (scaled=FALSE), center and scale using
    ## the median and the provided function 'scalefn'.
    if(!scaled) { ## Center and scale the data to (0, 1) - robustly
        x <- doScale(x, center=median, scale=scalefn)$x
    }

    nsets <- 6
    hsets <- matrix(integer(), h, nsets)

    ## Determine 6 initial estimates (ordering of obs)
    ## 1. Hyperbolic tangent of standardized data
    y1 <- tanh(x)
    R1 <- cor(y1)
    P <- eigen(R1, symmetric=TRUE)$vectors
    hsets[,1] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## 2. Spearmann correlation matrix
    R2 <- cor(x, method="spearman")
    P <- eigen(R2, symmetric=TRUE)$vectors
    hsets[,2] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## 3. Tukey normal scores
    y3 <- qnorm((apply(x, 2L, rank) - 1/3)/(n + 1/3))
    R3 <- cor(y3, use = "complete.obs")
    P <- eigen(R3, symmetric=TRUE)$vectors
    hsets[,3] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## 4. Spatial sign covariance matrix
    znorm <- sqrt(rowSums(x^2))
    ii <- znorm > .Machine$double.eps
    x.nrmd <- x
    x.nrmd[ii,] <- x[ii, ] / znorm[ii]
    SCM <- crossprod(x.nrmd)# / (n-1) not needed for e.vectors
    P <- eigen(SCM, symmetric=TRUE)$vectors
    hsets[,4] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## 5. BACON
    ind5 <- order(znorm)
    half <- ceiling(n/2)
    Hinit <- ind5[1:half]
    covx <- cov(x[Hinit, , drop=FALSE])
    P <- eigen(covx, symmetric=TRUE)$vectors
    hsets[,5] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## 6. Raw OGK estimate for scatter
    P <- ogkscatter(x, scalefn, only.P = TRUE)
    hsets[,6] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## Now combine the six pack :
    if(full.h) hsetsN <- matrix(integer(), n, nsets)
    for(k in 1:nsets) ## sort each of the h-subsets in *increasing* Mah.distances
    {
        xk <- x[hsets[,k], , drop=FALSE]
	svd <- classPC(xk, signflip=FALSE) # [P,T,L,r,centerX,meanvct] = classSVD(xk)
        if(svd$rank < p) ## FIXME: " return("exactfit")  "
            stop('More than half of the observations lie on a hyperplane.')
        score <- (x - rep(svd$center, each=n)) %*% svd$loadings
        ord <- order(mahalanobisD(score, FALSE, sqrt(abs(svd$eigenvalues))))
        if(full.h)
            hsetsN[,k] <- ord
        else hsets[,k] <- ord[1:h]
    }
    ## return
    if(full.h) hsetsN else hsets
} ## {r6pack}
