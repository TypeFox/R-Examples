####========== Pairwise methods for covariance / correlation =================


### From: Kjell Konis <konis@stats.ox.ac.uk>
### To: R-SIG-Robust@stat.math.ethz.ch, Ricardo Maronna ...
### Cc: Rand Wilcox ...
### Subject: Re: [RsR] [R] M-estimator R function question
### Date: Mon, 5 Dec 2005 10:29:11 +0000

### Here is an implementation of the OGK estimator completely in R.  I
### haven't touched it for a while and I forget how thoroughly I tested
### it so use it with a bit of caution.

###    http://www.stats.ox.ac.uk/~konis/pairwise.q
###    --------------------------------------------

##-------------------------------------------------------------------------------
## Computes the orthogonalized pairwise covariance matrix estimate described in
## in Maronna and Zamar (2002).

## Use: pairwise(X, 2, gk.sigmamu, gk, hard.rejection) for the
##      Gnanadesikan-Kettenring estimate.

## Alternatively, supply your own functions.

## MM  replaced   sweep(X, 1, .., '*') which is inefficient!
## ==  used       crossprod()  where appropriate
##
## I don't like the names gk.sigmamu() and gk(),
##      "gk":= Gnanadesikan-Kettenring; particularly not for the Tau-estimator
##      which is not at all related to G.-K.
##  ---> replacements   s/gk.sigmamu/scaleTau2/
##                      s/gk/covGK/
## -- also in the line of the many other cov*() functions I've renamed
##                      s/pairwise/covOGK/

## NOTA BENE: Is *now* consistent, since MM made  scaleTau2()  consistent

### Documentation -----> ../man/covOGK.Rd
##  =============        ================

##' Compute the mahalanobis distances for *diagonal* var/cov matrix:
##' @param x  n x p numeric data matrix
##' @param center numeric p-vector (or length 1 - is recycled) or FALSE
##' @param sd numeric p-vector of "standard deviations"
##' @examples all.equal(mahalanobisD(x, FALSE,    sd),
##'                     mahalanobis (x, rep(0,p), diag(sd^2)))
mahalanobisD <- function(x, center, sd) {
    ## Compute the mahalanobis distances (for diagonal cov).
    if(!identical(center, FALSE))
        x <- sweep(x, 2L, center, check.margin=FALSE)
    rowSums(sweep(x, 2L, sd, '/', check.margin=FALSE)^2)
}


covOGK <- function(X, n.iter = 2,
		   sigmamu,
		   rcov = covGK, weight.fn = hard.rejection,
		   keep.data = FALSE, ...)
{
    stopifnot(n.iter >= 1)
    call <- match.call()
    X <- as.matrix(X)

    p <- ncol(X)
    if(p < 2) stop("'X' must have at least two columns")

    Z <- X # as we use 'X' for the (re)weighting
    U <- diag(p)
    A <- list()

    ## Iteration loop.
    for(iter in 1:n.iter) { ## only a few iterations

	## Compute the vector of standard deviations d and
	## the covariance matrix U.

	d <- apply(Z, 2L, sigmamu, ...)
	Z <- sweep(Z, 2L, d, '/', check.margin=FALSE)

	for(i in 2:p) { # only need lower triangle of U
	    for(j in 1:(i - 1))
		U[i, j] <- rcov(Z[ ,i], Z[ ,j], ...)
	}

	## Compute the eigenvectors of U and store them as columns of E:
	## eigen(U, symmetric) only needs left/lower triangle
	E <- eigen(U, symmetric = TRUE)$vectors

	## Compute A and store it for each iteration
	A[[iter]] <- d * E

	## Project the data onto the eigenvectors
	Z <- Z %*% E
    }

    ## End of orthogonalization iterations.

    ## Compute the robust location and scale estimates for
    ## the transformed data.
    sqrt.gamma <- apply(Z, 2L, sigmamu, mu.too = TRUE, ...)
    center <- sqrt.gamma[1, ]
    sqrt.gamma <- sqrt.gamma[2, ]

    distances <- mahalanobisD(Z, center, sd=sqrt.gamma)

    ## From the inside out compute the robust location and
    ## covariance matrix estimates.  See equation (5).

    ## MM[FIXME]: 1st iteration (often the only one!) can be made *much* faster
    ##	  -----
    covmat <- diag(sqrt.gamma^2)

    for(iter in n.iter:1) {
	covmat <- A[[iter]] %*% covmat %*% t(A[[iter]])
	center <- A[[iter]] %*% center
    }

    center <- as.vector(center)

    ## Compute the reweighted estimate.	 First, compute the
    ## weights using the user specified weight function.

    weights <- weight.fn(distances, p, ...)
    sweights <- sum(weights)

    ## Then compute the weighted location and covariance
    ## matrix estimates.

    ## MM FIXME 2 : Don't need any of this, if all weights == 1
    ##	  -----	   (which is not uncommon) ==> detect that "fast"

    wcenter <- colSums(X * weights) / sweights
    Z <- sweep(X, 2L, wcenter, check.margin=FALSE) * sqrt(weights)
    wcovmat <- crossprod(Z) / sweights

    list(center = center,
	 cov = covmat,
	 wcenter = wcenter,
	 wcov = wcovmat,
	 weights = weights,
	 distances = distances,
	 n.iter = n.iter,
	 sigmamu = deparse(substitute(sigmamu)),
	 weight.fn = deparse(substitute(weight.fn)),
	 rcov = deparse(substitute(rcov)),
	 call = call,
	 ## data.name = data.name,
	 data = if(keep.data) X)
}


## a version with weights and consistency (but only one tuning const!!)
## is in /u/maechler/R/other-people/Mspline/Mspline/R/scaleTau.R
##
scaleTau2 <- function(x, c1 = 4.5, c2 = 3.0, consistency = TRUE,
                      mu.too = FALSE, ...)
{
    ## NOTA BENE: This is *NOT* consistency corrected
    n <- length(x)
    medx <- median(x)
    x. <- abs(x - medx)
    sigma0 <- median(x.) ## = MAD(x)  {without consistency factor}
    mu <-
        if(c1 > 0) {
            ## w <- pmax(0, 1 - (x. / (sigma0 * c1))^2)^2   -- but faster:
            x. <- x. / (sigma0 * c1)
            w <- 1 - x.*x.
            w <- ((abs(w) + w)/2)^2

            sum(x * w) / sum(w) }
        else medx

    x <- (x - mu) / sigma0
    rho <- x^2
    rho[rho > c2^2] <- c2^2
    ## sigma2 <- sigma0^2 * sum(rho)/ n

    if(!identical(consistency,FALSE)) {
	Erho <- function(b)
	    ## E [ rho_b ( X ) ]   X ~ N(0,1)
	    2*((1-b^2)*pnorm(b) - b * dnorm(b) + b^2) - 1
	Es2 <- function(c2)
	    ## k^2 * E[ rho_{c2} (X' / k) ] , where X' ~ N(0,1), k= qnorm(3/4)
	    Erho(c2 * qnorm(3/4))
        ## the asymptotic E[ sigma^2(X) ]  is Es2(c2):
        ## TODO: 'n-2' below will probably change; therefore not yet documented
        nEs2 <- (if(consistency == "finiteSample") n-2 else n) * Es2(c2)
    } else nEs2 <- n

    ## return
    c(if(mu.too) mu,
      ## sqrt(sigma2) == sqrt( sigma0^2 / n * sum(rho) ) :
      sigma0 * sqrt(sum(rho)/nEs2))
}

## Two other simple 'scalefun' to be used for covOGK;
## s_Qn(), s_Sn() are in ./qnsn.R
s_mad <- function(x, mu.too= FALSE, na.rm = FALSE) {
    if (na.rm) x <- x[!is.na(x)]
    mx <- median(x)
    c(if(mu.too) mx, mad(x, center = mx))
}
s_IQR <- function(x, mu.too= FALSE, na.rm = FALSE) {
    Qx <- quantile(x, (1:3)/4, na.rm = na.rm, names = FALSE)
    c(if(mu.too) Qx[2], (Qx[3] - Qx[1]) * 0.5 * formals(mad)$constant)
}

covGK <- function(x, y, scalefn = scaleTau2, ...)
{
    ## Gnanadesikan-Kettenring's, based on   4*Cov(X,Y) = Var(X+Y) - Var(X-Y)
    (scalefn(x + y, ...)^2 - scalefn(x - y, ...)^2) / 4
}

hard.rejection <- function(distances, p, beta = 0.9, ...)
{
    d0 <- median(distances) * qchisq(beta, p) / qchisq(0.5, p)
    wts <- double(length(distances))# == 0, but
    wts[distances <= d0] <- 1.0
    wts
}

##-- TODO "pairwise QC" ... etc
##--> ~maechler/R/MM/STATISTICS/robust/pairwise-new.R
