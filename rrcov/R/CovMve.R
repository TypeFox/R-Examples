CovMve <- function(x,
                   alpha=1/2,
                   nsamp=500,
                   seed=NULL,
                   trace=FALSE,
                   control)
{

    use.correction <- FALSE
    ## Analize and validate the input parameters ...

    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.

    defcontrol <- CovControlMcd()       # default control
    if(!missing(control)){
        if(alpha == defcontrol@alpha)       alpha <- control@alpha
        if(nsamp == defcontrol@nsamp)       nsamp <- control@nsamp
        if(is.null(seed) || seed == defcontrol@seed)         seed <- control@seed
        if(trace == defcontrol@trace)       trace <- control@trace
    }
    tolSolve <- defcontrol@tolSolve

    xcall <- match.call()

    if(length(seed) > 0) {
        if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
            seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
        }
        assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    ##   Options "best" and "exact" for nsamp
    ##   nsamp will be further analized in the wrapper .fastmcd()
    if(!missing(nsamp) && is.numeric(nsamp) && nsamp <= 0)
        stop("Invalid number of trials nsamp = ",nsamp, "!")

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
              dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
        stop("All observations have missing values!")
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]
    h <- h.alpha.n(alpha, n, p) # h(alpha) , the size of the subsamples
    ######## h <- floor(n/2)

    if(n <= p + 1)          # ==> floor((n+p+1)/2) > n - 1  -- not Ok
        stop(if (n <= p)    # absolute barrier!
             "n <= p -- you can't be serious!"
        else "n == p+1  is too small sample size for MCD")
    ## else
    if(n < 2 * p) { ## p+1 < n < 2p
        warning("n < 2 * p, i.e., possibly too small sample size")
    }
    if(h > n)
        stop("Sample size n  <  h(alpha; n,p) := size of \"good\" subsample")
    else if(alpha > 1)
        stop("alpha must be <= 1")

    ## VT::29.07.2008 - raw.cnp2 and cnp2 are vectors of size 2 and  will
    ##   contain the correction factors (concistency and finite sample)
    ##   for the raw and reweighted estimates respectively. Set them
    ##   initially to 1.  If use.correction is set to FALSE
    ##   (default=FALSE), the finite sample correction factor will not
    ##   be used (neither for the raw estimates nor for the reweighted)
    ##
    ##  The finite sample correction factors for MVE are not known, except
    ##  for the very old cnp2=(1 + 15/(n - p))^2 (Rousseeuw&van Zomeren, 1988)
    ##  therefore will remain always 1.
    ##
    ##  The consistancy factor for the raw covariance is
    ##      cf <- median(dist)/qchisq(0.5, p)
    ##
    ## FIXME: should be
    ##      cf <- quantile(dist, alpha)/qchisq(alpha, p)
    ##

    raw.cnp2 <- cnp2 <- c(1,1)

    ##  Case: alpha==1
    ##  ...
    ##  ...
    ##  Case p == 1

    method <- "Minimum volume ellipsoid estimator"
    mve <- .fastmve(x, h, nsamp)

    ## Compute the consistency correction factor for the raw MCD
    ##  (see calfa in Croux and Haesbroeck)
    ## calpha <- .MCDcons(p, h/n)    ## VT::19.3.2007
    mvecov <- cov.wt(x[mve$best,])
    rcenter <- mvecov$center
    rcov <- mvecov$cov
    mah <- mahalanobis(x, rcenter, rcov, tol = tolSolve)
    calpha <- quantile(mah, h/n)/qchisq(h/n, p)   # as in MASS
    names(calpha) <- NULL
    correct <- if(use.correction) (1 + 15/(n - p))^2 else 1.
    raw.cnp2 <- c(calpha, correct)

    rcov <- calpha * correct * rcov

    ## Again consider p == 1
    ## else, i.e. p >= 2
        ## handle exact fit, i.e. not general position situtations
        ##
        ## OK, in general position and mve$cov is not singular
        ##  do reweighting

        ## FIXME: here we assume that mve$cov is not singular
        ## ----- but it could be!
        quantiel <- qchisq(0.975, p)
        mah <- mahalanobis(x, rcenter, rcov, tol = tolSolve)
        weights <- as.numeric(mah < quantiel)
        sum.w <- sum(weights)

        ## Compute and apply the consistency correction factor for
        ## the reweighted cov

        if(sum.w == n) {
            cdelta.rew <- 1
            correct.rew <- 1
        }else {
            cdelta.rew <- .MCDcons(p, sum.w/n) ## VT::: 19.07.2008
            correct.rew <- if(use.correction) 1 else 1.
            cnp2 <- c(cdelta.rew, correct.rew)
        }

        xcov <- cov.wt(x, wt = weights)
        xcov$cov <- cdelta.rew * correct.rew * xcov$cov

        raw.mah <- mah
        raw.weights <- weights

        ## Check if the reweighted scatter matrix is singular and
        ##  compute distances and weights based on it
        if( - (determinant(xcov$cov, logarithm = TRUE)$modulus[1] - 0)/p > 50) {
            ## ans$singularity <- list(kind = "reweighted.MCD")
            if(trace)
                cat("The reweighted MCD scatter matrix is singular.\n")
            mah <- raw.mah
        }
        else {
            mah <- mahalanobis(x, xcov$center, xcov$cov, tol = tolSolve)
            weights <- as.numeric(mah < quantiel)
        }

    ##
    ans <- new("CovMve",
                call=xcall,
                iter=nsamp,
                crit=mve$scale,
                cov=xcov$cov,
                center=xcov$center,
                mah=mah,
                wt=weights,
                n.obs=n,
                X=x,
                method=method,
                best=mve$best,
                alpha=alpha,
                quan=h,
                raw.center=rcenter,
                raw.cov=rcov,
                raw.mah=raw.mah,
                raw.wt=raw.weights,
                raw.cnp2=raw.cnp2,
                cnp2=cnp2
            )
    ans
}

.fastmve <- function(x, h, nsamp)
{
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
    nind <- p+1

    tmp <- .C('r_fast_mve',
        x = if(is.double(x)) x else as.double(x),
        as.integer(n),
        as.integer(p),
        as.integer(nsamp),
        nsing = as.integer(0),
        ctr = as.double(rep(0,p)),
        cov = as.double(rep(0,p*p)),
        scale = as.double(0),
        best=as.integer(rep(0,n)),
        as.integer(nind),
        as.integer(h),
        as.double(qchisq(0.5, p)),
        PACKAGE="rrcov")

    mve.cov <- matrix(tmp$cov, p, p)
    return(list(center= tmp$ctr,
                cov=mve.cov,
                scale=tmp$scale,
                best=tmp$best[1:h],
                nsamp=nsamp,
                nsing = tmp$nsing))
}
