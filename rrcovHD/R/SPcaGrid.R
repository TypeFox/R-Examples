setMethod("getQuan", "SPcaGrid", function(obj) obj@n.obs)

##  The S3 version
SPcaGrid <- function (x, ...)
    UseMethod("SPcaGrid")

SPcaGrid.formula <- function (formula, data = NULL, subset, na.action, ...)
{
    cl <- match.call()

    mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0)
        stop("response not allowed in formula")
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    ## this is not a 'standard' model-fitting function,
    ## so no need to consider contrasts or levels
    if (.check_vars_numeric(mf))
        stop("PCA applies only to numerical variables")

    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0
    x <- model.matrix(mt, mf)

    res <- SPcaGrid.default(x, ...)

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("SPcaGrid")
    res@call <- cl

#    if (!is.null(na.act)) {
#        res$na.action <- na.act
#        if (!is.null(sc <- res$x))
#            res$x <- napredict(na.act, sc)
#    }

    res
}

SPcaGrid.default <- function(x, k=0, kmax=ncol(x),  method = c ("mad", "sd", "qn", "Qn"),
    lambda = 1, scale=FALSE, na.action = na.fail, trace=FALSE, ...)
{

    cl <- match.call()

    if(missing(x)){
        stop("You have to provide at least some data")
    }
    data <- as.matrix(x)
    n <- nrow(data)
    p <- ncol(data)

    ##
    ## verify and set the input parameters: k and kmax
    ##
    kmax <- max(min(floor(kmax), rankMM(x)),1)
    if(trace)
        cat("k=", k, ", kmax=", kmax, ".\n", sep="")

    if((k <- floor(k)) < 0)
        k <- 0
    else if(k > kmax) {
        warning(paste("The number of principal components k = ", k, " is larger than kmax = ", kmax, "; k is set to ", kmax,".", sep=""))
        k <- kmax
    }
    if(k != 0)
        k <- min(k, ncol(data))
    else
    {
        k <- min(kmax, ncol(data))
        if(trace)
            cat("The number of principal components is defined by the algorithm. It is set to ", k,".\n", sep="")
    }


    ######################################################################
    if(is.logical(scale))
    {
        scale <- if(scale) sd else  NULL
    }
    method <- match.arg(method)
    if(method== "Qn")
        method <- "qn"
    out <- sPCAgrid(x, k, lambda=lambda, method=method, scale=scale, ...)

    scores <- predict(out)
    center   <- out$center
    scale <- out$scale
    sdev     <- out$sdev
    scores   <- as.matrix(scores[, 1:k])
    loadings <- as.matrix(out$loadings[, 1:k])
    eigenvalues  <- (sdev^2)[1:k]

    ######################################################################
    names(eigenvalues) <- NULL
    if(is.list(dimnames(data)))
        rownames(scores) <- rownames(data)          # dimnames(scores)[[1]] <- dimnames(data)[[1]]
    dimnames(scores)[[2]] <- paste("PC", seq_len(ncol(scores)), sep = "")
    dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("SPcaGrid")
    res <- new("SPcaGrid", call=cl,
                            loadings=loadings,
                            eigenvalues=eigenvalues,
                            center=center,
                            scale=scale,
                            scores=scores,
                            k=k,
                            n.obs=n)

    ## Compute distances and flags
    res <- rrcov::pca.distances(res, x, p)
    return(res)
}

###############
##  utilities and help functions
##
##

##  calculate the degree of sparsity for the loadings vector ll
.dos <- function(ll, zero.tol=1e-2)
{

    q <- ncol(ll)
    ret <- vector(mode="numeric", length=q)
    for(i in 1:q)
        ret[i] <- length(which(abs(ll[,i]) < zero.tol))
    ret
}

## Calculate the CPEV of x for k PCs
##  CPEV=Cumulative Percentage of Explained Variance
.CPEV <- function(x, k=1)
{
    vars1 <- getEigenvalues(x);  vars1 <- vars1/sum(vars1)
    cvars1 <- cumsum(vars1)
    cvars1[k]
}

###
##  tradeoff() - extract sparse PC for different values of lambda and
##  calculate the corresponding explained variance
.tradeoff <- function(x, k, lambda.max=2.5, lambda.n=10, method="sd")
{
    p <- ncol(x)
    lambda.seq <- seq(0, lambda.max, length.out=lambda.n)
    cpev <- vector(mode="numeric", length=length(lambda.seq))
    i <- 1
    for(lambda in lambda.seq)
    {
        spc <- SPcaGrid(x, lambda=lambda, method=method, k=p)
        cpev[i] <- .CPEV(spc, k)
        cat("\n", i, round(lambda, 2), cpev[i], .dos(getLoadings(spc))[1:k], "\n")
        i <- i+1
    }

    ret <- cbind.data.frame(lambda.seq, cpev)
    colnames(ret) <- c("lambda", "CPEV")
    ret
}

.check_vars_numeric <- function(mf)
{
    ## we need to test just the columns which are actually used.
    mt <- attr(mf, "terms")
    mterms <- attr(mt, "factors")
    mterms <- rownames(mterms)[apply(mterms, 1, any)]
    any(sapply(mterms, function(x) is.factor(mf[,x]) || !is.numeric(mf[,x])))
}
