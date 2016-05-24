setMethod("getCenter", "Pca", function(obj) obj@center)
setMethod("getScale", "Pca", function(obj) obj@scale)
setMethod("getLoadings", "Pca", function(obj) obj@loadings)
setMethod("getEigenvalues", "Pca", function(obj) obj@eigenvalues)
setMethod("getSdev", "Pca", function(obj) sqrt(obj@eigenvalues))
setMethod("getScores", "Pca", function(obj) obj@scores)
setMethod("getPrcomp", "Pca", function(obj) {
    ret <- list(sdev=sqrt(obj@eigenvalues),
         rotation=obj@loadings,
         center=obj@center,
         scale=obj@scale,
         x=obj@scores)
    class(ret) <- "prcomp"
    ret
})

##
## Follow the standard methods: show, print, plot
##
setMethod("show", "Pca", function(object) myPcaPrint(object))
setMethod("summary", "Pca", function(object, ...){
    vars <- getEigenvalues(object)
    vars <- vars/sum(vars)
    importance <- rbind("Standard deviation" = getSdev(object),
                        "Proportion of Variance" = round(vars,5),
                        "Cumulative Proportion" = round(cumsum(vars), 5))
    colnames(importance) <- colnames(getLoadings(object))
    new("SummaryPca", pcaobj=object, importance=importance)
})
setMethod("show", "SummaryPca", function(object){

    cat("\nCall:\n")
    print(object@pcaobj@call)

    digits = max(3, getOption("digits") - 3)

    cat("Importance of components:\n")
    print(object@importance, digits = digits)
    invisible(object)
})

## setMethod("print",     "Pca", function(x, ...) myPcaPrint(x, ...))

setMethod("predict",   "Pca", function(object, ...){
    predict(getPrcomp(object), ...)
})
setMethod("screeplot", "Pca", function(x, ...){
    screeplot(getPrcomp(x), ...)
})

setMethod("biplot",    "Pca", function(x, choices=1L:2L, scale=1, ...){

    if(length(getEigenvalues(x)) < 2)
        stop("Need at least two components for biplot.")

    lam <- sqrt(getEigenvalues(x)[choices])
    scores <- getScores(x)
    n <- NROW(scores)
    lam <- lam * sqrt(n)

    if(scale < 0 || scale > 1)
        warning("'scale' is outside [0, 1]")

    if(scale != 0)
        lam <- lam^scale
    else
        lam <- 1

    xx <- t(t(scores[, choices]) / lam)
    yy <- t(t(getLoadings(x)[, choices]) * lam)
    .biplot(xx, yy, ...)
    invisible()
})

setMethod("scorePlot", "Pca", function(x, i=1, j=2, ...){
    pca.scoreplot(obj=x, i=i, j=j, ...)
})

##  The __outlier map__ (diagnostic plot, distance-distance plot)
##  visualizes the observations by plotting their orthogonal
##  distance to the robust PCA subspace versus their robust
##  distances within the PCA subspace. This allows to classify
##  the data points  into 4 types: regular observations, good
##  leverage points, bad leverage points and orthogonal outliers.
##  The outlier plot is only possible when k < r (the number of
##  selected components is less than the rank of the matrix).
##  Otherwise a __distance plot__ will be shown (distances against
##  index).
##
##  The __screeplot__ shows the eigenvalues and is helpful to select
##  the number of principal components.
##
##  The __biplot__ is plot which aims to represent both the
##  observations and variables of a matrix of multivariate data
##  on the same plot.
##
## The __scoreplot__ shows a scatterplot of i-th against j-th score
##  of the Pca object with superimposed tollerance (0.975) ellipse
##
## VT::17.06.2008
##setMethod("plot", "Pca", function(x, y="missing",
setMethod("plot", signature(x="Pca", y="missing"), function(x, y="missing",
                                id.n.sd=3,
                                id.n.od=3,
                                ...){

    if(all(x@od > 1.E-06))
        pca.ddplot(x, id.n.sd, id.n.od, ...)
    else
        pca.distplot(x, id.n.sd, ...)
})

myPcaPrint <- function(x, print.x=FALSE, print.loadings=FALSE, ...) {
    if(!is.null(cl <- x@call)) {
        cat("Call:\n")
        dput(cl)
        cat("\n")
    }

    cat("Standard deviations:\n"); print(sqrt(getEigenvalues(x)), ...)
    if(print.loadings)
    {
        cat("\nLoadings:\n");          print(getLoadings(x), ...)
    }

    if (print.x) {
        cat("\nRotated variables:\n"); print(getScores(x), ...)
    }
    invisible(x)
}

## Internal function to calculate the score and orthogonal distances and the
##  appropriate cutoff values for identifying outlying observations
##
##  obj  - the Pca object
##  data -
##  r    - rank
##  crit - criterion for computing cutoff for SD and OD
##
##  - cutoff for score distances: sqrt(qchisq(crit, k)
##  - cutoff for orthogonal distances: Box (1954)
##
pca.distances <- function(obj, data, r, crit=0.975) {
    .distances(data, r, obj, crit)
}
.distances <- function(data, r, obj, crit=0.975) {

    ## remember the criterion, could be changed by the user
    obj@crit.pca.distances <- crit

    ## compute the score distances and the corresponding cutoff value
    n <- nrow(data)
    smat <- diag(obj@eigenvalues, ncol=ncol(obj@scores))

    ## VT::02.06.2010: it can happen that the rank of the matrix
    ##  is nk=ncol(scores), but the rank of the diagonal matrix of
    ##  eigenvalues is lower: for example if the last singular
    ##  value was 1E-7, the last eigenvalue will be sv^2=1E-14
    ##
    nk <- min(ncol(obj@scores), rankMM(smat))
    if(nk < ncol(obj@scores))
        warning(paste("Too small eigenvalue(s): ", obj@eigenvalues[ncol(obj@scores)], "- the diagonal matrix of the eigenvalues cannot be inverted!"))

    obj@sd <- sqrt(mahalanobis(as.matrix(obj@scores[,1:nk]), rep(0, nk), diag(obj@eigenvalues[1:nk], ncol=nk)))
    obj@cutoff.sd <- sqrt(qchisq(crit, obj@k))

    ## Compute the orthogonal distances and the corresponding cutoff value
    ##  For each point this is the norm of the difference between the
    ##  centered data and the back-transformed scores
##    obj@od <- apply(data - repmat(obj@center, n, 1) - obj@scores %*% t(obj@loadings), 1, vecnorm)
    obj@od <- apply(data - matrix(rep(obj@center, times=n), nrow=n, byrow=TRUE) - obj@scores %*% t(obj@loadings), 1, vecnorm)

    if(is.list(dimnames(obj@scores))) {
        names(obj@od) <- dimnames(obj@scores)[[1]]
    }

    ## The orthogonal distances make sence only if the number of PCs is less than
    ##  the rank of the data matrix - otherwise set it to 0
    obj@cutoff.od <- 0
    if(obj@k != r) {
        obj@cutoff.od <- .crit.od(obj@od, crit=crit, classic=inherits(obj,"PcaClassic"))
    }

    ## flag the observations with 1/0 if the distances are less or equal the
    ##  corresponding  cutoff values
    obj@flag <- obj@sd <= obj@cutoff.sd
    if(obj@cutoff.od > 0)
        obj@flag <- (obj@flag & obj@od <= obj@cutoff.od)

    return(obj)
}

.crit.od <- function(od, crit=0.975, umcd=FALSE, quan, classic=FALSE)
{
    od <- od^(2/3)
    if(classic)
    {
        t <- mean(od)
        s <- sd(od)
    }else if(umcd)
    {
        ms <- unimcd(od, quan=quan)
        t <- ms$tmcd
        s <- ms$smcd
    }else
    {
        t <- median(od)
        s <- mad(od)
    }
    cv <- (t + s * qnorm(crit))^(3/2)
    cv
}

## Flip the signs of the loadings
##  - comment from Stephan Milborrow
##
.signflip <- function(loadings)
{
    if(!is.matrix(loadings))
        loadings <- as.matrix(loadings)
    apply(loadings, 2, function(x) if(x[which.max(abs(x))] < 0) -x else x)
}

##' @title Compute Classical Principal Components via SVD or Eigen
##' @param x a numeric matrix
##' @param scale logical
##' @param signflip
##' @param via.svd logical indicating if the SVE
##' @return a list
##' @author Valentin Todorov; efficiency tweaks by Martin Maechler
classPC <- function(x, scale=FALSE, signflip=TRUE, via.svd = n > p)
{
    if(!is.numeric(x) || !is.matrix(x))
        stop("'x' must be a numeric matrix")
    else if((n <- nrow(x)) <= 1)
        stop("The sample size must be greater than 1 for svd")
    p <- ncol(x)
    center <- colMeans(x)
    x <- scale(x, center=center, scale=scale)
    ##   -----
    if(isTRUE(scale))
        scale <- attr(x, "scaled:scale")

    if(via.svd) {
        svd <- svd(x/sqrt(n-1), nu=0)
        rank <- rankMM(x, sv=svd$d)
        loadings <- svd$v[,1:rank]
        eigenvalues <- (svd$d[1:rank])^2 ## FIXME: here .^2; later sqrt(.)
    } else { ## n <= p; was "kernelEVD"
        e <- eigen(tcrossprod(x)/(n-1), symmetric=TRUE)
        tolerance <- n * max(e$values) * .Machine$double.eps
        rank <- sum(e$values > tolerance)
        ii <- seq_len(rank)
        eigenvalues <- e$values[ii]
        ## MM{FIXME (efficiency)}:
        loadings <- t((x/sqrt(n-1))) %*% e$vectors[,1:rank] %*% diag(1/sqrt(eigenvalues))
    }

    ## VT::15.06.2010 - signflip: flip the sign of the loadings
    if(signflip)
        loadings <- .signflip(loadings)

    list(loadings=loadings,
         ## scores = x %*% loadings,
         eigenvalues=eigenvalues,
         rank=rank, center=center, scale=scale)
}

## VT::15.06.2010 - Added scaling and flipping of the loadings
##
classSVD <- function(x, scale=FALSE, signflip=TRUE){
    if(!is.numeric(x) || !is.matrix(x))
        stop("'x' must be a numeric matrix")
    else if(nrow(x) <= 1)
        stop("The sample size must be greater than 1 for svd")

    n <- nrow(x)
    p <- ncol(x)

    center <- apply(x, 2, mean)
    x <- scale(x, center=TRUE, scale=scale)
    if(scale)
        scale <- attr(x, "scaled:scale")

    svd <- svd(x/sqrt(n-1))

    rank <- rankMM(x, sv=svd$d)
    eigenvalues <- (svd$d[1:rank])^2
    loadings <- svd$v[,1:rank]

    ## VT::15.06.2010 - signflip: flip the sign of the loadings
    if(!is.matrix(loadings))
        loadings <- data.matrix(loadings)
    if(signflip)
        loadings <- .signflip(loadings)
    scores <- x %*% loadings

    list(loadings=loadings,
         scores=scores,
         eigenvalues=eigenvalues,
         rank=rank,
         center=center,
         scale=scale)
}

## VT::15.06.2010 - Added scaling and flipping of the loadings
##
kernelEVD <- function(x, scale=FALSE, signflip=TRUE){
    if(!is.numeric(x) || !is.matrix(x))
        stop("'x' must be a numeric matrix")
    else if(nrow(x) <= 1)
        stop("The sample size must be greater than 1 for svd")

    n <- nrow(x)
    p <- ncol(x)

    if(n > p) classSVD(x, scale=scale, signflip=signflip)
    else {
        center <- apply(x, 2, mean)
        x <- scale(x, center=TRUE, scale=scale)
        if(scale)
            scale <- attr(x, "scaled:scale")

        e <- eigen(x %*% t(x)/(n-1))

        tolerance <- n * max(e$values) * .Machine$double.eps
        rank <- sum(e$values > tolerance)

        eigenvalues <- e$values[1:rank]
        loadings <- t((x/sqrt(n-1))) %*% e$vectors[,1:rank] %*% diag(1/sqrt(eigenvalues))

        ## VT::15.06.2010 - signflip: flip the sign of the loadings
        if(signflip)
            loadings <- .signflip(loadings)

        scores <- x %*% loadings
        ret <- list(loadings=loadings,
                    scores=scores,
                    eigenvalues=eigenvalues,
                    rank=rank,
                    center=center,
                    scale=scale)
    }
}

## Score plot of the Pca object 'obj' - scatterplot of ith against jth score
##  with superimposed tollerance (0.975) ellipse
pca.scoreplot <- function(obj, i=1, j=2, main, id.n=0, ...)
{
    if(missing(main))
    {
        main <- if(inherits(obj,"PcaClassic")) "Classical PCA" else "Robust PCA"
    }

    x <- cbind(getScores(obj)[,i], getScores(obj)[,j])

## VT::11.06.2012
##  Here we assumed that the scores are not correlated and
##  used a diagonal matrix with the eigenvalues on the diagonal to draw the ellipse
##  This is not the case with PP methods, therefore we compute the covariance of
##  the scores, considering only the non-outliers
##      (based on sore and orthogonal distances)
##
##    ev <- c(getEigenvalues(obj)[i], getEigenvalues(obj)[j])
##    cpc <- list(center=c(0,0), cov=diag(ev), n.obs=obj@n.obs)

    flag <- obj@flag
    cpc <- cov.wt(x, wt=flag)

## We need to inflate the covariance matrix with the proper size:
##  see Maronna et al. (2006), 6.3.2, page 186
##
##  - multiply the covariance by quantile(di, alpha)/qchisq(alpha, 2)
##  where alpha = h/n
##
##    mdx <- mahalanobis(x, cpc$center, cpc$cov)
##    alpha <- length(flag[which(flag != 0)])/length(flag)
##    cx <- quantile(mdx, probs=alpha)/qchisq(p=alpha, df=2)
##    cpc$cov <- cx * cpc$cov


    .myellipse(x, xcov=cpc,
        xlab=paste("PC",i,sep=""),
        ylab=paste("PC",j, sep=""),
        main=main,
        id.n=id.n, ...)
    abline(v=0)
    abline(h=0)
}

## Distance-distance plot (or diagnostic plot, or outlier map)
## Plots score distances against orthogonal distances
pca.ddplot <- function(obj, id.n.sd=3, id.n.od=3, main, xlim, ylim, off=0.02, ...) {

    if(missing(main))
    {
        main <- if(inherits(obj,"PcaClassic")) "Classical PCA" else "Robust PCA"
    }


    if(all(obj@od <= 1.E-06))
        warning("PCA diagnostic plot is not defined")
    else
    {
        if(missing(xlim))
            xlim <- c(0, max(max(obj@sd), obj@cutoff.sd))
        if(missing(ylim))
            ylim <- c(0, max(max(obj@od), obj@cutoff.od))

        plot(obj@sd, obj@od, xlab="Score distance",
                             ylab="Orthogonal distance",
                             main=main,
                             xlim=xlim, ylim=ylim, type="p", ...)
        abline(v=obj@cutoff.sd)
        abline(h=obj@cutoff.od)
        label.dd(obj@sd, obj@od, id.n.sd, id.n.od, off=off)
    }
    invisible(obj)
}

## Distance plot, plots score distances against index
pca.distplot <- function(obj, id.n=3, title, off=0.02, ...) {

    if(missing(title))
    {
        title <- if(inherits(obj,"PcaClassic")) "Classical PCA" else "Robust PCA"
    }

    ymax <- max(max(obj@sd), obj@cutoff.sd)
    plot(obj@sd, xlab="Index", ylab="Score distance", ylim=c(0,ymax), type="p", ...)
    abline(h=obj@cutoff.sd)
    label(1:length(obj@sd), obj@sd, id.n, off=off)
    title(title)
    invisible(obj)
}

label <- function(x, y, id.n=3, off=0.02){
    xrange <- par("usr")
    xrange <- xrange[2] - xrange[1]
    if(id.n > 0) {
        n <- length(y)
        ind <- sort(y, index.return=TRUE)$ix
        ind <- ind[(n-id.n+1):n]
        if(is.character(names(y)))
            lab <- names(y[ind])
        else
            lab <- ind
        text(x[ind] - off*xrange, y[ind], lab)
    }
}

label.dd <- function(x, y, id.n.sd=3, id.n.od=3, off=0.02){
    xrange <- par("usr")
    xrange <- xrange[2] - xrange[1]
    if(id.n.sd > 0 && id.n.od > 0) {
        n <- length(x)
        ind.sd <- sort(x, index.return=TRUE)$ix
        ind.sd <- ind.sd[(n - id.n.sd + 1):n]
        ind.od <- sort(y, index.return=TRUE)$ix
        ind.od <- ind.od[(n - id.n.od + 1):n]
        lab <- ind.od
        if(is.character(names(y)))
            lab <- names(y[ind.od])
        text(x[ind.od] - off*xrange, y[ind.od], lab)
        lab <- ind.sd
        if(is.character(names(x)))
            lab <- names(x[ind.sd])
        text(x[ind.sd] - off*xrange, y[ind.sd], lab)
    }
}

## VT::30.09.2009 - add a parameter 'classic' to generate a default caption
##  "Robust biplot" or "Classical biplot" for a robust/classical
##  PCA object, resp.
##  --- do not use it for now ---
##
.biplot <- function(x, y, classic,
                    var.axes = TRUE,
                    col,
                    cex = rep(par("cex"), 2),
                    xlabs = NULL, ylabs = NULL,
                    expand=1,
                    xlim = NULL, ylim = NULL,
                    arrow.len = 0.1,
                    main = NULL, sub = NULL, xlab = NULL, ylab = NULL, ...)
{
    n <- nrow(x)
    p <- nrow(y)

##    if(is.null(main))
##        main <- if(classic) "Classical biplot" else "Robust biplot"

    if(missing(xlabs))
    {
       xlabs <- dimnames(x)[[1L]]
       if(is.null(xlabs))
            xlabs <- 1L:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
    if(missing(ylabs))
    {
       ylabs <- dimnames(y)[[1L]]
       if(is.null(ylabs))
            ylabs <- paste("Var", 1L:p)
    }
    ylabs <- as.character(ylabs)
    dimnames(y) <- list(ylabs, dimnames(y)[[2L]])

    if(length(cex) == 1L)
        cex <- c(cex, cex)

    pcol <- par("col")
    if(missing(col))
    {
       col <- par("col")
       if(!is.numeric(col))
            col <- match(col, palette(), nomatch=1L)
       col <- c(col, col + 1L)
    }
    else if(length(col) == 1L)
        col <- c(col, col)

    unsigned.range <- function(x)
        c(-abs(min(x, na.rm=TRUE)), abs(max(x, na.rm=TRUE)))
    rangx1 <- unsigned.range(x[, 1L])
    rangx2 <- unsigned.range(x[, 2L])
    rangy1 <- unsigned.range(y[, 1L])
    rangy2 <- unsigned.range(y[, 2L])

    if(missing(xlim) && missing(ylim))
    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if(missing(xlim)) xlim <- rangx1
    else if(missing(ylim)) ylim <- rangx2
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    on.exit(par(op))
    op <- par(pty = "s")
    if(!is.null(main))
        op <- c(op, par(mar = par("mar")+c(0,0,1,0)))
    plot(x, type = "n", xlim = xlim, ylim = ylim, col = col[[1L]],
         xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
    text(x, xlabs, cex = cex[1L], col = col[[1L]], ...)
    par(new = TRUE)
    plot(y, axes = FALSE, type = "n", xlim = xlim*ratio, ylim = ylim*ratio,
     xlab = "", ylab = "", col = col[[1L]], ...)

##    axis(3, col = col[2L], ...)
##    axis(4, col = col[2L], ...)
##    box(col = col[1L])
    axis(3, col = pcol, ...)
    axis(4, col = pcol, ...)
    box(col = pcol)

    text(y, labels=ylabs, cex = cex[2L], col = col[[2L]], ...)
    if(var.axes)
    arrows(0, 0, y[,1L] * 0.8, y[,2L] * 0.8, col = col[[2L]], length=arrow.len)
    invisible()
}
