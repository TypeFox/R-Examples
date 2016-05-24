setMethod("getQuan", "PcaNA", function(obj) obj@n.obs)


.getCov <- function(pc, k=ncol(getLoadings(pc)), inverse=FALSE)
{
    P <- getLoadings(pc)
    k <- min(ncol(P), k)

    lambda <- if(inverse) diag(1/getSdev(pc)[1:k]^2) else diag(getSdev(pc)[1:k]^2)

    P[, 1:k] %*% lambda %*% t(P[, 1:k])
}

##  The S3 version
PcaNA <- function (x, ...)
    UseMethod("PcaNA")

PcaNA.formula <- function (formula, data = NULL, subset, na.action, ...)
{
    cl <- match.call()

    mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0)
        stop("response not allowed in formula")
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    ## this is not a `standard' model-fitting function,
    ## so no need to consider contrasts or levels
    if (.check_vars_numeric(mf))
        stop("PCA applies only to numerical variables")

    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0
    x <- model.matrix(mt, mf)

    res <- PcaNA.default(x, ...)

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("PcaNA")
    res@call <- cl

#    if (!is.null(na.act)) {
#        res$na.action <- na.act
#        if (!is.null(sc <- res$x))
#            res$x <- napredict(na.act, sc)
#    }

    res
}

PcaNA.default <- function(x, k=ncol(x), kmax=ncol(x), conv=1e-10, maxiter=100,
    method=c("cov", "locantore", "hubert", "grid", "proj", "class"), cov.control=NULL,
    scale=FALSE, signflip=TRUE, crit.pca.distances=0.975, trace=FALSE, ...)
{
##    conv = 1e-10;
##    maxiter = 1000;

    cl <- match.call()
    method <- match.arg(method)

    if(missing(x))
    {
        stop("You have to provide at least some data")
    }

    data <- as.matrix(x)
    n <- nrow(data)
    p <- ncol(data)

    isna <- is.na(x)
    complete <- length(which(isna)) == 0

    xc <- .scaleNA(x)           # center and scale
    xc$x <- .inimiss(xc$x)
    critold <- sum(xc$x[isna]^2)

    for(iter in 1:maxiter)
    {
        pca <- if(method == "grid" | method == "proj") .PcaRobust(xc$x, k=k, control=method, ...) else .PcaRobust(xc$x, k=k, control=method, cov.control=cov.control, ...)
        xp <- getScores(pca)[, 1:k] %*% t(getLoadings(pca)[,1:k])
        xc$x[isna] <- xp[isna]      # replace missing elements with
                                    # new estimates
        crit <- sum(xc$x[isna]^2)
        if(trace)
            cat("\niter, crit, critold: ", iter, crit, critold, crit-critold, abs(crit - critold)/critold, "\n")

        if(critold == 0 || abs(crit - critold)/critold <= conv)
            break

        critold = crit
    }

    x <- .scaleNABack(xc$x, loc=xc$loc, sc=xc$sc)

    pca <- .PcaRobust(x, k=k, control=method, ...)

######################################################################
#    names(eigenvalues) <- NULL
#    if(is.list(dimnames(data)))
#    {
#        rownames(scores) <- rownames(data)  # dimnames(scores)[[1]] <- dimnames(data)[[1]]
#    }
#    dimnames(scores)[[2]] <- paste("PC", seq_len(ncol(scores)), sep = "")
#   dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))

    ## fix up call to refer to the generic, but leave arg name as 'formula'
#    cl[[1]] <- as.name("PcaCov")

    res <- new("PcaNA", call=cl,
                            loadings=pca@loadings,
                            eigenvalues=pca@eigenvalues,
                            center=pca@center,
                            scale=pca@scale,
                            scores=pca@scores,
                            k=pca@k,
                            n.obs=pca@n.obs,
                            Ximp=as.matrix(x))

    ## Compute distances and flags
    res <- pca.distances(res, x, p)

    return(res)
}

##
##  Center and scale the data matrix x
##
.scaleNA <- function(x)
{
    ## center and scale
    loc <- colMeans(x, na.rm = TRUE)

##  VT::deprecated sd(<data.frame>)
##    sc  <- sd(x, na.rm = TRUE)
  sc  <- apply(x, 2, sd, na.rm = TRUE)

    x <- sweep(x, 2, loc, "-")
    x <- sweep(x, 2, sc, "/")

    list(x=x, loc=loc, sc=sc)
}

##
##  Transform back a data matrix x with given location and scale
##
.scaleNABack <- function(x, loc, sc)
{
    ## back transform
    x <- sweep(x, 2, sc, "*")
    x <- sweep(x, 2, loc, "+")
    x
}

##  Initial estimation of missing elements
##
##      By default the missing elements of X are estimated as the
##      average of the corresponding row and column means.
##      If the optional parameter 'const' is supplied, the missing
##      elements are replaced by 'const'
##
##  INPUT:
##      X   - data matrix X (objects x variables) with missing elements
##
##  OUTPUT:
##      X   - data matrix with initial estimates of missing elements
##
.inimiss <- function(x, const)
{
    isna <- is.na(x)
    if(length(which(isna)) == 0)
        return (x)

    if(!missing(const))
    {
        x[isna] <- const
        return(x)
    }

    n <- nrow(x)
    p <- ncol(x)
    cmn <- colMeans(x, na.rm=TRUE)
    rmn <- rowMeans(x, na.rm=TRUE)

    xx <- matrix(0,nrow=n,ncol=p)
    for(i in 1:n)
        for(j in 1:p)
            xx[i,j] <- (rmn[i] + cmn[j])/2
   isna <- is.na(x)
   x[isna] <- xx[isna]
   x
}

##  Move this function to rrcov - similarly to CovRobust
##
##  control can be a character specifying the name of the estimate, one of:
##  locantore, hubert, proj, grid, cov, class
##  If no control object is given or 'auto' is selected, the choice of the
##  estimator will depend on the size of the data:
##  - to be defined
##
.PcaRobust <- function(x, control, ...)
{
    cl <- match.call()
    if(missing(x)){
        stop("You have to provide at least some data")
    }
    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))

    n <- nrow(x)
    p <- ncol(x)

    mm <- NULL
    if(missing(control))
    {
        mm <- "auto"
    }else if(is.character(control))
    {
        mm <- casefold(control)
    }

    ## either no control specified or the estimator is given by a character name -
    ##  create the necessary control object.
    if(!is.null(mm))
    {
        pca <- switch(mm,
            auto = {
                PcaCov(x, ...)
            },
            locantore = PcaLocantore(x, ...),
            hubert = PcaHubert(x, ...),
            proj = PcaProj(x, ...),
            grid = PcaProj(x, ...),
            cov  = PcaCov(x, ...),
            class = PcaClassic(x, ...))

        ## this is the 'default' option of the switch
        if(is.null(pca))
            stop(paste("Undefined estimator: ", mm, "- Must be one of [locantore, hubert, proj, grid, cov, class]"))
    } else
            stop(paste("Undefined estimator: ", mm, "- Must be one of [locantore, hubert, proj, grid, cov, class]"))
    pca
}

.check_vars_numeric <- function(mf)
{
    ## we need to test just the columns which are actually used.
    mt <- attr(mf, "terms")
    mterms <- attr(mt, "factors")
    mterms <- rownames(mterms)[apply(mterms, 1, any)]
    any(sapply(mterms, function(x) is.factor(mf[,x]) || !is.numeric(mf[,x])))
}
