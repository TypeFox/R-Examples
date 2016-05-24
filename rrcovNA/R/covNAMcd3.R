.covNAMcd <- function(x,
           cor = FALSE,
           alpha = 1/2,
           nsamp = 500,
           seed = NULL,
           trace = FALSE,
           use.correction = TRUE,
           imputation = TRUE,
           impMeth = c("norm" , "seq", "rseq"),
           control = rrcov.control())
{
    ##   Analyze and validate the input parameters ...

    ## if a control object was supplied, take the option parameters
    ## from it, but if single parameters were passed (not defaults)
    ## they will override the control object.
    ## defCtrl <- rrcov.control()# default control
    # 'control' specified
    if(missing(alpha))  alpha <- control$alpha
    if(missing(nsamp))  nsamp <- control$nsamp
    if(missing(seed))    seed <- control$seed
    if(missing(trace))  trace <- control$trace
    if(missing(use.correction)) use.correction <- control$use.correction

    tolSolve <- control$tolSolve # had 1e-10 hardwired {now defaults to 1e-14}

    if(length(seed) > 0) {
        if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
            seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
        }
        assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    ##   vt::03.02.2006 - added options "best" and "exact" for nsamp
    ##   nsamp will be further analized in the wrapper .fastmcd()
    if(!missing(nsamp) && is.numeric(nsamp) && nsamp <= 0)
        stop("Invalid number of trials nsamp = ",nsamp, "!")

    impMeth <- match.arg(impMeth)

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
                    dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows which contain only missings
    na.x <- rowSums(ifelse(is.na(x),1,0)) == ncol(x)
    ok <- !na.x
    x <- x[ok, , drop = FALSE]

    dx <- dim(x)
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]

    ## h(alpha) , the size of the subsamples
    h <- h.alpha.n(alpha, n, p)
    if(n <= p + 1) # ==> floor((n+p+1)/2) > n - 1  -- not Ok
        stop(if (n <= p) "n <= p -- you can't be serious!"
             else "n == p+1  is too small sample size for MCD")
    ## else
    if(n < 2 * p) { ## p+1 < n < 2p
        warning("n < 2 * p, i.e., possibly too small sample size")
        ## was stop("Need at least 2*(number of variables) observations ")
    }

    if(h > n)
        stop("Sample size n  <  h(alpha; n,p) := size of \"good\" subsample")
    else if(alpha > 1)
        stop("alpha must be <= 1")

    cutoff <- 0.975
    quantiel <- qchisq(cutoff, p)
    raw.cnp2 <- cnp2 <- c(1,1)
    ans <- list(method = "Minimum Covariance Determinant Estimator for incomplete data.", call = match.call())

## CASE 1: alpha=1 - Just compute the classical estimates --------
    if(alpha == 1) {
        ## mcd <- cov.wt(x)$cov
        ## loc <- as.vector(colMeans(x))
        cc <- .cov.na.wt(x)
        mcd <- cc$cov
        loc <- cc$center
        obj <- determinant(mcd, logarithm = TRUE)$modulus[1]
        if ( -obj/p > 50 ) {
            ans$cov <- mcd
            dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
            if (cor)
                ans$cor <- cov2cor(ans$cov)
            ans$center <- loc
            if(length(dimn[[2]]))
                names(ans$center) <- dimn[[2]]
            ans$n.obs <- n
            ans$singularity <- list(kind = "classical")
            if(trace) cat("classical estimate is singular\n")

            weights <- 1
        }
        else {
##            mah <- mahalanobis(x, loc, mcd, tol = tolSolve)
##            weights <- as.numeric(mah < quantiel) # 0/1
##            sum.w <- sum(weights)
##            ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
            xcov <- .cov.na.rew(x, loc, mcd, tol = tolSolve)
            weights <- xcov$wt
            sum.w <- sum(weights)

            ans$center <- xcov$center
            ans$cov <- xcov$cov
            ans$n.obs <- xcov$n.obs
######################            ans$cov <- sum.w/(sum.w - 1) * ans$cov

            ## Consistency factor for reweighted MCD
            if(sum.w != n) {
                cnp2[1] <- .MCDcons(p, sum.w/n)
                ans$cov <- ans$cov * cnp2[1]
            }
            if( - (determinant(ans$cov, logarithm = TRUE)$modulus[1] - 0)/p > 50) {
                ans$singularity <- list(kind = "reweighted.MCD")
                if(trace) cat("reweighted MCD is singular\n")
            }
            else {
##                mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
##                weights <- as.numeric(mah < quantiel) # 0/1
                mah <- .mah.na(x, ans$center, ans$cov, tol = tolSolve)
                if(!is.null(mah))
                {
                    ans$mah <- mah$d
                    weights <- as.integer(.dflag(mah, pr=cutoff))
                }else
                {
                    ans$mah <- ans$raw.mah
                }

            }
        }

        ans$alpha <- alpha
        ans$quan <- h
        ans$raw.cov <- mcd
        ans$raw.center <- loc
        if(!is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$crit <- exp(obj)
        ans$method <-
            paste(ans$method,
                  "\nThe minimum covariance determinant estimates based on",
                  n, "observations \nare equal to the classical estimates.")

        ans$mcd.wt <- weights
        if(length(dimn[[1]]))
            names(ans$mcd.wt) <- dimn[[1]]
        ans$wt <- NULL
        ans$X <- x
        if(length(dimn[[1]]))
            dimnames(ans$X)[[1]] <- names(ans$mcd.wt)
        else
            dimnames(ans$X) <- list(seq(1,n), NULL)
        if(trace)
            cat(ans$method, "\n")
        ans$raw.cnp2 <- raw.cnp2
        ans$cnp2 <- cnp2
        class(ans) <- "mcd"
        return(ans)
    } ## end {alpha=1} --

## CASE 2: alpha<1 - Compute MCD --------

    if(p == 1) {
        stop("Univariate estimation with missing data is not supported!")
    } ## end p=1

    if(imputation){
        mcd <- .cov.na.imp(x, alpha=alpha, nsamp=nsamp, seed=seed, impMeth)
    }
    else {
        mcd <- .fastmcdNA(x, h, nsamp)
        if(mcd$exactfit == 3)
            mcd <- .cov.na.imp(x, alpha=alpha, nsamp=nsamp, seed=seed, impMeth)
    }

    ## Compute the consistency correction factor for the raw MCD
    ##  (see calfa in Croux and Haesbroeck)
    calpha <- .MCDcons(p, h/n)    ## VT::19.3.2007
    correct <- if(use.correction) .MCDcnp2(p, n, alpha) else 1.
    raw.cnp2 <- c(calpha, correct)

      ## Apply correction factor to the raw estimates
      ## and use them to compute weights
      mcd$initcovariance <- calpha * mcd$initcovariance * correct
      dim(mcd$initcovariance) <- c(p, p)

      if(mcd$exactfit != 0) {
        ans$cov <- mcd$initcovariance
        ans$center <- as.vector(mcd$initmean)
        if(!is.null(nms <- dimn[[2]])) {
            dimnames(ans$cov) <- list(nms, nms)
            names(ans$center) <- nms
        }
        ans$n.obs <- n

        if(!(mcd$exactfit %in% c(1,2,3)))
            stop("Unexpected 'exactfit' code ", mcd$exactfit, ". Please report!")

        ans$singularity <- list(kind = "on.hyperplane", exactCode = mcd$exactfit)
        if(mcd$exactfit == 3){
            ans$singularity = "Missing data"
        }

        ans$alpha <- alpha
        ans$quan <- h
        ans$raw.cov <- mcd$initcovariance
        ans$raw.center <- as.vector(mcd$initmean)
        if(!is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        weights <- rep(0,n)
        ans$crit <- 0
      } ## end exact fit <==>  (mcd$exactfit != 0)
      else { ## exactfit == 0 : have general position ------------------------

##        mah <- mahalanobis(x, mcd$initmean, mcd$initcovariance, tol = tolSolve)
##        mcd$weights <- weights <- as.numeric(mah < quantiel)

        ## Compute the reweighted matrix taking into account that the data are
        ##  incomplete. .cov.na.rew() returns also the weights with which the reweighting
        ##  was done
        xcov <- .cov.na.rew(x, mcd$initmean, mcd$initcovariance, tol = tolSolve)
        weights <- xcov$wt
        sum.w <- sum(weights)

        ## Compute and apply the consistency correction factor for
        ## the reweighted cov
        if(sum.w == n) {
            cdelta.rew <- 1
            correct.rew <- 1
        }
        else {
            cdelta.rew <- .MCDcons(p, sum.w/n) ## VT::19.3.2007
            correct.rew <- if(use.correction) .MCDcnp2.rew(p, n, alpha) else 1.
            cnp2 <- c(cdelta.rew, correct.rew)
        }

        ans$center <- xcov$center
        ans$cov <- xcov$cov
        ans$n.obs <- xcov$n.obs

#        ans <- c(ans, cov.wt(x, wt = weights, cor))
#        ans$cov <- sum.w/(sum.w - 1) * ans$cov

        ans$cov <- ans$cov * cdelta.rew * correct.rew
        if(!is.null(nms <- dimn[[2]])) {
            dimnames(ans$cov) <- list(nms, nms)
            names(ans$center) <- nms
        }

        ans$best <- sort(as.vector(mcd$best))
        ans$alpha <- alpha
        ans$quan <- h
        ans$raw.cov <- mcd$initcovariance
        ans$raw.center <- as.vector(mcd$initmean)
        if(!is.null(nms <- dimn[[2]])) {
            names(ans$raw.center) <- nms
            dimnames(ans$raw.cov) <- list(nms,nms)
        }
        ans$raw.weights <- weights
        ans$crit <- mcd$mcdestimate
##        ans$raw.mah <- mahalanobis(x, ans$raw.center, ans$raw.cov, tol = tolSolve)
        mah <- .mah.na(x, ans$raw.center, ans$raw.cov, tol = tolSolve)
        ans$raw.mah <- mah$d

        ## Check if the reweighted scatter matrix is singular.
        if( - (determinant(ans$cov, logarithm = TRUE)$modulus[1] - 0)/p > 50) {
            ans$singularity <- list(kind = "reweighted.MCD")
            if(trace) cat("The reweighted MCD scatter matrix is singular.\n")
            ans$mah <- ans$raw.mah
        }
        else {
            ## mah <- mahalanobis(x, ans$center, ans$cov, tol = tolSolve)
            ## ans$mah <- mah
            ## weights <- as.numeric(mah < quantiel)
            mah <- .mah.na(x, ans$center, ans$cov, tol = tolSolve)
            if(!is.null(mah))
            {
                ans$mah <- mah$d
                weights <- as.integer(.dflag(mah, pr=cutoff))
            }else
            {
                ans$mah <- ans$raw.mah
            }
        }
      } ## end{ not exact fit }

    ans$mcd.wt <- weights
    if(length(dimn[[1]]))
        names(ans$mcd.wt) <- dimn[[1]]
    ans$wt <- NULL
    ans$X <- x
    if(length(dimn[[1]]))
        dimnames(ans$X)[[1]] <- names(ans$mcd.wt)
    else
        dimnames(ans$X) <- list(seq(1,n), NULL)
    if(trace)
        cat(ans$method, "\n")
    ans$raw.cnp2 <- raw.cnp2
    ans$cnp2 <- cnp2
    class(ans) <- "mcd"
    return(ans)
}

.mymedian <- function(x) median(x, na.rm=TRUE)
.mymad <- function(x) mad(x, na.rm=TRUE)

.fastmcdNA <- function(x, h, nsamp){

    rescale <- TRUE

    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

if(rescale){
    xmed <- apply(x, 2, .mymedian)       # claculate MEDIAN for each variable
    xmad <- apply(x, 2, .mymad)          # claculate MADN for each variable
    x <- sweep(x, 2, xmed, "-")
    x <- sweep(x, 2, xmad, "/")
}
    xmcd <- .fastmcdNAx(x, h, nsamp)

if(rescale){
    xmcd$initcovariance <- matrix(xmcd$initcovariance, nrow=p)
    xmcd$initmean <- xmcd$initmean * xmad + xmed
    xmcd$initcovariance <- t(t(xmad * xmcd$initcovariance) * xmad)
    xmcd$mcdestimate <- xmcd$mcdestimate * prod(xmad * xmad)
    xmcd$initcovariance <- as.vector(xmcd$initcovariance)
}


## an ML estimate of cov is returned - convert to an unbiased
    nsel <- length(xmcd$best)
    fc <- nsel/(nsel-1)
    xmcd$initcovariance <- xmcd$initcovariance * fc
    xmcd$mcdestimate <- xmcd$mcdestimate * fc ^ p
##    print(fc)
##    print(matrix(xmcd$initcovariance,ncol=p))
##    print(xmcd$mcdestimate)
    xmcd
}

.fastmcdNAx <- function(x, h, nsamp)
{

    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

    ##   parameters for partitioning {equal to those in Fortran !!}
    kmini <- 5
    nmini <- 300
    km10 <- 10*kmini
    nmaxi <- nmini*kmini

    ##   vt::03.02.2006 - added options "best" and "exact" for nsamp
    if(!missing(nsamp)) {
    if(is.numeric(nsamp) && (nsamp < 0 || nsamp == 0 && p > 1)) {
        warning("Invalid number of trials nsamp= ", nsamp,
                    " ! Using default.\n")
        nsamp <- -1
    } else if(nsamp == "exact" || nsamp == "best") {
        myk <- p + 1 ## was 'p'; but p+1 ("nsel = nvar+1") is correct
        if(n > 2*nmini-1) {
                ## TODO: make 'nmini' user configurable; however that
                ##      currently needs changes to the fortran source
                ##  and re-compilation!
        warning("Options 'best' and 'exact' not allowed for n greater than ",
            2*nmini-1,".\nUsing default.\n")
        nsamp <- -1
        } else {
        nall <- choose(n, myk)
        if(nall > 5000 && nsamp == "best") {
            nsamp <- 5000
            warning("'nsamp = \"best\"' allows maximally 5000 subsets;\n",
                "computing these subsets of size ",
                            myk," out of ",n,"\n")
        } else { ## "exact" or ("best"  &  nall < 5000)
            nsamp <- 0 ## all subsamples
            if(nall > 5000)
            warning("Computing all ", nall, " subsets of size ",myk,
                " out of ",n,
                "\n This may take a very long time!\n",
                immediate. = TRUE)
        }
        }
    }

    if(!is.numeric(nsamp) || nsamp == -1) { # still not defined - set it to the default
        defCtrl <- rrcov.control() # default control
        if(!is.numeric(nsamp))
        warning("Invalid number of trials nsamp= ",nsamp,
            " ! Using default nsamp= ",defCtrl$nsamp,"\n")
        nsamp <- defCtrl$nsamp  # take the default nsamp
    }
    }

    ##  Dimensions of work arrays for handling the mising values
    ##  xint(nint), xdble(ndble) and xdat(n*p)
    d <- (2 + 3 * p + p^2)/2
    nint <- (p+1)*(p+1) + n*p + 3*n + 3*p
    ndble <- 4*d + 3*p + n*p
    mvcode <-  -99999

    x[is.na(x)] <- as.double(mvcode)

    ##   Allocate temporary storage for the Fortran implementation,
    ##   directly in the .Fortran() call.
    ##    (if we used C, we'd rather allocate there, and be quite faster!)

    .Fortran("FNANMCD",
         x = if(is.double(x)) x else as.double(x),
         n =    as.integer(n),
         p =    as.integer(p),   ## = 'nvar'  in Fortran
         nhalff =   as.integer(h),
         nsamp  =   as.integer(nsamp),# = 'krep'
         initcovariance = double(p * p),
         initmean       = double(p),
         best       = rep.int(as.integer(10000), h),
         mcdestimate = double(1),    ## = 'det'
         weights   = integer(n),
         exactfit  = integer(1), # output indicator: 0: ok; 1: ..., 2: ..
         coeff     = matrix(double(5 * p), nrow = 5, ncol = p), ## plane
         kount     = integer(1),
         adjustcov = double(p * p), ## << never used -- FIXME
         integer(1),## << 'seed' no longer used -- FIXME
         temp   = integer(n),
         index1 = integer(n),
         index2 = integer(n),
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

         cstock = double(10 * p * p),# (10,nvmax2)
         mstock = double(10 * p),    # (10,nvmax)
         c1stock = double(km10 * p * p), # (km10,nvmax2)
         m1stock = double(km10 * p), # (km10,nvmax)
         dath = double(nmaxi * p),   # (nmaxi,nvmax)

         cutoff = qchisq(0.975, p),
         chimed = qchisq(0.5,   p),
         xdat = double(n * p),
         d = as.integer(d),
         xint = integer(nint),
         nint = as.integer(nint),
         xdble = double(ndble),
         ndble = as.integer(ndble),
         mvcode = as.double(mvcode),
         xindex1 = integer(n),
         xindex2 = integer(n), PACKAGE="rrcovNA")[
         c("initcovariance", "initmean", "best", "mcdestimate", "exactfit") ]
}

.cov.na.imp <- function(x, alpha = 1/2, nsamp = 500, seed = NULL, impMeth){

    if(impMeth == "norm")
    {
        ## impute the missing data using package norm
        s <- prelim.norm(x)                     # do preliminary manipulations
        thetahat <- em.norm(s, showits=FALSE)   # find the mle
        rngseed(1234567)                        # set random number generator seed
        ximp <- imp.norm(s, thetahat, x)        # impute missing data under the MLE
        xx<-imp.norm(s, thetahat, x)            # impute missing data under the MLE
    }else if(impMeth == "seq")
    {
        ximp <- impSeq(x)
    }else if(impMeth == "rseq")
    {
        ximp <- impSeqRob(x, alpha=0.75)$x
    }else
        stop("Undefined imputation method in .cov.na.imp(): ", impMeth)

    mcd <- covMcd(ximp, alpha=alpha, nsamp=nsamp, seed=seed, trace=FALSE)

    exactfit <- if(!is.null(mcd$singularity)) 2 else 0

##  the raw covariance was adjusted - restore it

    mcd$raw.cov <- mcd$raw.cov / prod(mcd$raw.cnp2)

    list(initcovariance=mcd$raw.cov,
               initmean=mcd$raw.center,
                   best=mcd$best,
            mcdestimate=mcd$crit,
               exactfit=exactfit)

}
