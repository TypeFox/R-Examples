.xmiss <- function(x)
{
    n <- nrow(x)
    p <- ncol(x)
    ret <- vector(n, mode="logical")
    for(i in 1:n)
    {
        ret[i] <- length(which(!is.na(x[i,]))) < p
    }
    ret
}

.cov.na.wt <- function(x, wt = rep(1, nrow(x)), method = c("unbiased", "ML"))
{

    method <- match.arg(method)
    if(is.data.frame(x))
        x <- as.matrix(x)
    else if(!is.matrix(x))
        stop("'x' must be a matrix or a data frame")

    n <- nrow(x)
    if(with.wt <- !missing(wt)) {
        if(length(wt) != n)
            stop("length of 'wt' must equal the number of rows in 'x'")
        if(any(wt < 0) || (s <- sum(wt)) == 0)
            stop("weights must be non-negative and not all zero")
    }

    cc <- .cov.na(x[which(wt == 1),])

    ## this is an ML estimate - scale the cov if "unbiased" requested
    sum.w <- sum(wt)
    if(method == "unbiased")
        cc$cov <- sum.w/(sum.w - 1) * cc$cov
    cc$wt <- wt
    cc
}

.cov.na <- function(x){

    mvcode <- -99999

## get dimensions of x
    n <- nrow(x)
    p <- ncol(x)
    storage.mode(x) <- "double"
    x <- .na.to.snglcode(x, mvcode)

    d <- (2 + 3 * p + p^2)/2
    nint <- (p+1)*(p+1) + n*p + 3*n + 3*p
    ndble <- 4*d + 3*p + n*p

    mu <- double(p)
    sigma <- matrix(0, p, p)

    start <- .Fortran("emnint",
                x = if(is.double(x)) x else as.double(x),
                as.integer(n), as.integer(p), as.integer(d),
                integer(nint), as.integer(nint),
                double(ndble), as.integer(ndble),
                center=mu,
                cov=sigma,
                as.double(mvcode),
                PACKAGE="rrcovNA")[c("center", "cov")]

    names(start$center) <- dimnames(x)[[2]]
    dimnames(start$cov) <- list(names(start$center), names(start$center))
    start$n.obs <- n
    start
}

.cov.na.rew <- function(x, center, cov, cutoff=0.975, method = c("unbiased", "ML"), ...){

    method <- match.arg(method)
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
    mah <- .mah.na(x, center, cov, ...)
    weights <- .dflag(mah, pr=cutoff)

    .cov.na.wt(x, as.integer(weights), method)
}

##
## Flag observations as outliers based on the mahalanobis distances 'd'
##  using qchisq(pr, p).
##
##  If p is not supplied, it is expected that the mahalanobis distances
##  are computed from incomplete data and 'd' contains a vector 'pp'
##  with the number of non-missing values in each observation.
##
.dflag <- function(d, pr=0.975, p){
    if(is.null(d))
        return(NULL)

    if(missing(p) && !is.null(d$pp))
        retval <- sqrt(d$d) < sqrt(qchisq(pr, d$pp))
    else if(!missing(p))
    {
        chi <- qchisq(pr, p)
        retval <- sqrt(d) < sqrt(chi)
    } else
        stop("Error - invalid call")
    retval
}

##  Compute mahalanobis distance for incomplete data
##  Returns a list of two vectors:
##  d - squared mahalanobis distances and
##  pp - number of observed values in each observation
##      to be used as degrees of freedom in chisq
##
.mah.na <- function (x, center, cov, inverted=FALSE, ...)
{
    x <- if (is.vector(x))
            matrix(x, ncol = length(x))
         else
            as.matrix(x)
    x <- sweep(x, 2, center)

    n <- nrow(x)
    p <- ncol(x)
    d <- pp <- z <- vector(mode="numeric", length=n)
    cinv1 <- cov
    if (!inverted) {
        cinv1 <- try(solve(cov, ...))
        if(!is.numeric(cinv1))
        {
          print(cov)
          cat("\n============ Singular cov matrix in .mah.na Returning NULL =======\n")
          return (NULL)
          ## stop(.Last.value)
        }
    }
    for(i in 1:n)
    {
        indo <- which(!is.na(x[i,]))
        pp[i] <- length(indo)
        if(pp[i] == p)   # complete observation - use cinv1
        {
            d[i] <- rowSums((x[i,indo] %*% cinv1) * x[i,indo])
        } else
        {
            cinv <- try(solve(cov[indo,indo], ...))
            if(!is.numeric(cinv))
            {
            print(indo)
            print(cov)
            print(cov[indo,indo])
            }
            d[i] <- rowSums((x[i,indo] %*% cinv) * x[i,indo])
        }
        z[i] <- .wh(d[i], pp[i])
    }
    names(d) <- rownames(x)
    names(pp) <- rownames(x)
    names(z) <- rownames(x)
    list(d=d, pp=pp, z=z)
}

## Wilson-Hilferty transformation
.wh <- function(d,pp)
    ((d/pp)^(1/3) - 1 + 2/(9*pp))/sqrt(2/(9*pp))
