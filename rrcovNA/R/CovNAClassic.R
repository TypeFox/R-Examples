## Just a shortcut and for competability with older versions
##  same as CovNAClassic()
##
CovNA <- function(x, unbiased = TRUE)
{
    return (CovNAClassic(x, unbiased))
}

CovNAClassic <- function(x, unbiased = TRUE)
{

    if(is.data.frame(x))
    {
        x <- data.matrix(x)
    } else if(!is.matrix(x))
    {
        x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))
    }

    call <- match.call()
    method = "EM algorithm for incomplete normal data."

    ## drop all rows which contain only missings
    na.x <- rowSums(ifelse(is.na(x),1,0)) == ncol(x)
    ok <- !na.x
    x <- x[ok, , drop = FALSE]

    dimn <- dimnames(x)
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

    s <- prelim.norm(x)                     # do preliminary manipulations
    thetahat <- em.norm(s, showits=FALSE)   # find the mle estimates
    xcov <- getparam.norm(s, thetahat)      # extract the estimated location and covariance
    if(!is.null(nms <- dimn[[2]]))
    {
        dimnames(xcov$sigma) <- list(nms, nms)
        names(xcov$mu) <- nms
    }


    ans <- new("CovNAClassic",
               call = call,
               cov=xcov$sigma,
               center=xcov$mu,
               n.obs=n,
               X = x,
               method=method)

    ans
}

## VT::17.06.2008
##setMethod("plot", "CovClassic", function(x, y="missing",
setMethod("plot", signature(x="CovNAClassic", y="missing"), function(x, y="missing",
                                which=c("all", "distance", "qqchi2", "tolEllipsePlot", "screeplot"),
                                ask = (which=="all" && dev.interactive(TRUE)),
                                cutoff,
                                id.n,
                                tol = 1e-7, ...)
{

    data <- getData(x)
    ##  parameters and preconditions
    if(is.vector(data) || is.matrix(data)) {
        if(!is.numeric(data))
            stop( "x is not a numeric dataframe or matrix.")
    } else if(is.data.frame(data)) {
        if(!all(sapply(data,data.class) == "numeric"))
            stop( "x is not a numeric dataframe or matrix.")
    }

    n <- dim(data)[1]
    p <- dim(data)[2]

    if(length(getCenter(x))  == 0 ||  length(getCov(x)) == 0)
        stop( "Invalid object: attributes center and cov missing!")

    if(length(getCenter(x))  != p)
        stop( "Data set and provided center have different dimensions!")

    ## Check for singularity of the cov matrix
    if(isSingular(x))
        stop("The covariance matrix is singular!")

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(!missing(id.n) && !is.null(id.n)) {
        id.n <- as.integer(id.n)
        if(id.n < 0 || id.n > n)
            stop(sQuote("id.n")," must be in {1,..,",n,"}")
    }

    mymiss <- .xmiss(x@X)

    md <- sqrt(getDistance(x)$d)
    which <- match.arg(which)

    op <- if (ask) par(ask = TRUE) else list()
    on.exit(par(op))

    ## index plot of mahalanobis distances
    if(which == "all" || which == "distance") {
        .mydistplotNA(md, mymiss, cutoff, classic=TRUE, id.n=id.n, ...)
    }

    ## qq-plot of the mahalanobis distances versus the
    ## quantiles of the chi-squared distribution
    if(which == "all" || which == "qqchi2") {
        .qqplotNA(md, p, mymiss, cutoff=cutoff, classic=TRUE, id.n=id.n, ...)
    }

##    if(which == "all" || which == "tolEllipsePlot") {
##        if(p == 2)
##            .tolellipse(ccov = x, cutoff=cutoff, id.n=id.n, tol=tol, ...)
##        else if(which != "all")
##            warning("Warning: For tolerance ellipses the dimension must be 2!")
##    }

    if(which == "all" || which == "screeplot") {
        myscreeplot(ccov=x)
    }
}) ## end { plot("CovClassic") }
