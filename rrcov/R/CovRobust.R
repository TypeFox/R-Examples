##  control can be a character specifying the name of the estimate, one of:
##  auto, mcd, ogk, m, mve, sfast, surreal, bisquare, rocke
##  If no control object is given or 'auto' is selected, the choice of the
##  estimator will depend on the size of the data:
##  - Stahel-Donoho: n < 1000 and p < 10 or n < 5000 and p < 5
##  - MCD: n < 50000 and p < 20
##  - OGK: otherwise
##
CovRobust <- function(x, control, na.action = na.fail)
{
    x <- na.action(x)
    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))

    n <- nrow(x)
    p <- ncol(x)

    control <- .covRobustControl(control, n, p)
    restimate(control, x)
}

.covRobustControl <- function(method, n, p)
{
    mm <- NULL
    if(missing(method))
        mm <- "auto"
    else if(is.character(method))
        mm <- casefold(method)

    ## either no control specified or the estimator is given by a character name -
    ##  create the necessary control object.
    if(!is.null(mm)){
        control <- switch(mm,
            auto = {
                if(missing(n) || missing(p))
                    CovControlMcd()
                else if(n < 1000 && p < 10 || n < 5000 && p < 5)
                    CovControlSde()
                else if(n < 50000 && p < 10)
                    CovControlSest(method="bisquare")
                else if(n < 50000 && p < 20)
                    CovControlSest(method="rocke")
                else
                    CovControlOgk(smrob="s_mad", svrob="qc")
            },
            sde = CovControlSde(),
            mcd = CovControlMcd(),
            ogk = CovControlOgk(),
            m   = CovControlMest(),
            mve = CovControlMve(),
            sfast = CovControlSest(method="sfast"),
            surreal = CovControlSest(method="surreal"),
            bisquare = CovControlSest(method="bisquare"),
            rocke = CovControlSest(method="rocke"))

        ## this is the 'default' option of the switch
        if(is.null(control))
            stop(paste("Undefined estimator: ", method))
    } else
        control <- method
    control
}

setMethod("isClassic", "CovRobust", function(obj) FALSE)
setMethod("isClassic", "SummaryCovRobust", function(obj) FALSE)
setMethod("getMeth", "CovRobust", function(obj) obj@method)
setMethod("getRaw", "CovRobust", function(obj){

    if(is(obj, "CovMcd") | is(obj, "CovMve") | is(obj, "CovOgk"))
    {
        obj@center <- obj@raw.center
        obj@cov <- obj@raw.cov
        obj@mah <- obj@raw.mah
        obj@wt <- obj@raw.wt
    }
    if(is(obj, "CovMcd") | is(obj, "CovMve"))
    {
        obj@cnp2 <- obj@raw.cnp2
    }
    invisible(obj)
})

##
## Follow the standard methods: show, summary, plot
##
setMethod("show", "CovRobust", function(object){
    cat("\nCall:\n")
    print(object@call)
    cat("-> Method: ", object@method, "\n")
    if(is.list(object@singularity))
        cat(strwrap(.MCDsingularityMsg(object@singularity, object@n.obs)), sep ="\n")

    digits = max(3, getOption("digits") - 3)
    cat("\nRobust Estimate of Location: \n")
    print.default(format(getCenter(object), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nRobust Estimate of Covariance: \n")
    print.default(format(getCov(object), digits = digits), print.gap = 2, quote = FALSE)
    invisible(object)
})

setMethod("summary", "CovRobust", function(object, ...){

    new("SummaryCovRobust", covobj=object, evals=eigen(object@cov)$values)

})


setMethod("show", "SummaryCovRobust", function(object){

    cat("\nCall:\n")
    print(object@covobj@call)

    digits = max(3, getOption("digits") - 3)
    cat("\nRobust Estimate of Location: \n")
    print.default(format(getCenter(object), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nRobust Estimate of Covariance: \n")
    print.default(format(getCov(object), digits = digits), print.gap = 2, quote = FALSE)

    cat("\nEigenvalues of covariance matrix: \n")
    print.default(format(getEvals(object), digits = digits), print.gap = 2, quote = FALSE)

    cat("\nRobust Distances: \n")
    print.default(format(as.vector(getDistance(object)), digits = digits), print.gap = 2, quote = FALSE)
})

## VT::17.06.2008
##setMethod("plot", "CovRobust", function(x, y="missing",
setMethod("plot", signature(x="CovRobust", y="missing"),
                                function(x, y="missing",
                                    which=c("dd", "all", "distance", "qqchi2",
                                    "tolEllipsePlot", "pairs", "screeplot", "xydistance", "xyqqchi2"),
                                classic= FALSE,
                                ask = (which=="all" && dev.interactive(TRUE)),
                                cutoff,
                                id.n,
                                labels.id = rownames(x$X),
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

    ccov <- CovClassic(data)
    md <- rd <- NULL
    if(!isSingular(ccov))
        md <- sqrt(getDistance(ccov))
    if(!isSingular(x))
        rd <- sqrt(getDistance(x))

    which <- match.arg(which)
    op <- if (ask) par(ask = TRUE) else list()
    on.exit(par(op))

    ## distance-distance plot: here we need both robust and mahalanobis distances
    if((which == "all" || which == "dd") && !is.null(md) && !is.null(rd)) {
        .myddplot(md, rd, cutoff=cutoff, id.n=id.n, ...) # distance-distance plot
    }

    ## index plot of mahalanobis distances
    if((which == "all" || which == "distance")  && !is.null(rd)) {
        ylim <- NULL
        if(classic && !is.null(md)) {
            opr <- if(prod(par("mfrow")) == 1) par(mfrow=c(1,2), pty="m") else list()

            ##VT::10.11.2007 - set same scale on both plots
            ylim <- c(min(rd,md), max(md,rd))
        }

        .mydistplot(rd, cutoff=cutoff, id.n=id.n, ...)            # index plot of robust distances
        if(classic && !is.null(md)) {
            .mydistplot(md, cutoff=cutoff, classic=TRUE, id.n=id.n, ylim=ylim, ...) # index plot of mahalanobis distances
            par(opr)
        }
    }

    ## lattice: index plot of mahalanobis distances
    if(which == "xydistance"  && !is.null(rd)) {
        print(.xydistplot(x, cutoff=cutoff, ...))      # lattice:  index plot of robust distances

    }

    ## qq-plot of the mahalanobis distances versus the
    ## quantiles of the chi-squared distribution
    if((which == "all" || which == "qqchi2")  && !is.null(rd)) {
        if(classic && !is.null(md)) {
            opr <- if(prod(par("mfrow")) == 1) par(mfrow=c(1,2), pty="m") else list()
        }
        .qqplot(rd, p, cutoff=cutoff, id.n=id.n, ...) # qq-plot of the robust distances versus the
                                                      # quantiles of the chi-squared distribution
        if(classic && !is.null(md)) {
            .qqplot(md, p, cutoff=cutoff, classic=TRUE, id.n=id.n, ...)
                                                 # qq-plot of the mahalanobis distances
            par(opr)
        }
    }

    ## lattice: qq-plot of the mahalanobis distances versus the
    ##          quantiles of the chi-squared distribution
    if(which == "xyqqchi2"  && !is.null(rd)) {
        print(.xyqqchi2(x, cutoff=cutoff, ...))        # lattice:  qq-plot of the distances versus
    }

    if(which == "tolEllipsePlot" || which == "pairs") {
        if(which == "tolEllipsePlot" & length(dim(data)) >= 2 && dim(data)[2] == 2){
            if(!is.null(rd)){
                if(classic &&  !is.null(md))
                    .tolellipse(rcov=x, ccov = ccov, cutoff=cutoff, id.n=id.n, tol=tol, ...)
                else
                    .tolellipse(rcov=x, cutoff=cutoff, id.n=id.n, tol=tol, ...)
            }
        }else if(length(dim(data)) >= 2 && dim(data)[2] <= 10)
        {
            .rrpairs(x, ...)
        }else if(which != "all")
            warning("Warning: For tolerance ellipses the dimension must be less than 10!")
    }

    if(which == "all" || which == "screeplot") {
        myscreeplot(ccov=ccov, rcov=x)
    }
}) ## end { plot("CovRobust") }
