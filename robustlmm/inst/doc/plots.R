require(ggplot2)
require(Matrix)

## image() for sparse Matrices from the Matrix package
## using ggplot methods

matrix2data <- function(x) {
    x <- as(x, "matrix")
    data <- data.frame(x=as.vector(x), expand.grid(i=1:nrow(x), j=1:ncol(x)))
    ## set 0 entries to zero
    data[data$x == 0,"x"] <- NA
    ## set dimensions
    attr(data, "nrow") <- nrow(x)
    attr(data, "ncol") <- ncol(x)
    data
}

image.ggplot <- function(x, limits, low = "#132B43", high = "#56B1F7",
                         data = matrix2data(x), ...) {
    if (missing("limits"))
        limits <- c(min(data$x, 0, na.rm=TRUE), max(data$x, 2, na.rm=TRUE))
    nrow <- attr(data, "nrow")
    ncol <- attr(data, "ncol")
    # levelplot(x ~ j + i, data)
    ggplot(data, aes(j, i, fill = x)) + geom_tile(...) + coord_equal(ratio = 1) +
        scale_fill_continuous(limits=limits, low=low, high=high, na.value = NA) +
            scale_x_continuous(expand = c(0, 0), limits = c(0.5, ncol+0.5)) +
                scale_y_reverse(expand = c(0, 0), limits = c(nrow+0.5, -0.5)) +
                    theme(legend.position="none") +
                        ylab("Row") + xlab("Column")
  }

image.ggplot.fast <- function(x, ...) {
    x <- as(x, "dgTMatrix")
    data <- data.frame(x=x@x, j=x@j + 1L, i=x@i + 1L)
    attr(data, "nrow") <- nrow(x)
    attr(data, "ncol") <- ncol(x)
    image.ggplot(..., data=data)
}

g.get_colors_brewer <- function(n, name='Dark2') {
    idx <- 1:n
    if (name=='Dark2' & n == 8) {
        idx <- c(6:3,7,1,2,8)[idx]
    }
    RColorBrewer::brewer.pal(n, name)[idx]
}
## plot(1:8, rep(1, 8), col = col[c(6:3,7,1:2,8)], pch = 18, cex = 5)

augLines <- function(obj, ..., bindAlso, subset) {
    data <- cbind(obj@frame, fitted=fitted(obj))
    if (!missing(bindAlso)) data <- cbind(data, bindAlso)
    env <- parent.frame()
    if (!missing(subset)) data <- data[eval(substitute(subset), data, env),]
    geom_line(data=data, aes(y=fitted), ...)
}
Mmat <- function(object) { ## create Matrix M
    lM <- object@pp$M()
    p <- object@pp$p
    q <- object@pp$q
    M <- Matrix(0, p+q, p+q)
    M[1:p, 1:p] <- lM$M_BB
    M[1:p, p+(1:q)] <- t(lM$M_bB)
    M[p+(1:q), 1:p] <- lM$M_bB
    M[p+(1:q), p+(1:q)] <- lM$M_bb
    M
}
Dmat <- function(object) { ## create augmented design matrix
    with(object@pp, {
        D <- Matrix(0, p+q, n+q)
        D[1:p, 1:n] <- t(.U_eX)
        D[p+(1:q), 1:n] <- t(.U_eZ %*% U_b)
        D[p+(1:q), n+(1:q)] <- -Lambda_b
        D })
}
vcovFull <- function(object) {
    ## unscaled V
    lC <- with(object@pp, bdiag(Diagonal(x=rep(Epsi2_e, n)), Epsi_bpsi_bt))
    MD <- Mmat(object) %*% Dmat(object)
    Vunsc <- MD %*% tcrossprod(lC, MD)
    V <- getME(object, "sigma")^2 * Vunsc
}
## all.equal(as(vcovFull(rfm31)[1:4, 1:4], "matrix"), as(vcov(rfm31), "matrix"),
##           check.attr=FALSE) ## close enough??
augLinesSE <- function(obj, ..., alpha=0.2, bindAlso, subset) {
    if (!is(obj, "rlmerMod")) obj <- as(obj, "rlmerMod")
    ## TRUE: do not include random effects
    data <- cbind(obj@frame, fitted={
        if (FALSE) obj@pp$X %*% fixef(obj) else fitted(obj)
        })
    lD <- with(obj@pp, {
        lD <- Matrix(0, n, p+q)
        lD[1:n, 1:p] <- X
        lD[1:n, p+(1:q)] <- crossprod(Zt, U_b)
        lD })
    se <- sqrt(diag(lD %*% tcrossprod(vcovFull(obj), lD)))
    ## add se bands
    data$fittedLow <- data$fitted - 1.96 * se
    data$fittedUp <- data$fitted + 1.96 * se
    if (!missing(bindAlso)) data <- cbind(data, bindAlso)
    env <- parent.frame()
    if (!missing(subset)) data <- data[eval(substitute(subset), data, env),]
    geom_ribbon(data=data, aes(ymin=fittedLow, ymax=fittedUp), alpha=alpha, ...)
}
augLinesPop <- function(obj, ..., bindAlso, subset) {
    if (!is(obj, "rlmerMod")) obj <- as(obj, "rlmerMod")
    data <- cbind(obj@frame, fitted=obj@pp$X %*% fixef(obj))
    if (!missing(bindAlso)) data <- cbind(data, bindAlso)
    env <- parent.frame()
    if (!missing(subset)) data <- data[eval(substitute(subset), data, env),]
    geom_line(data=data, aes(y=fitted), ...)
}
augLinesPopSE <- function(obj, ..., alpha = 0.2, bindAlso, subset) {
    if (!is(obj, "rlmerMod")) obj <- as(obj, "rlmerMod")
    X <- obj@pp$X
    se <- sqrt(diag(tcrossprod(X %*% vcov(obj), X)))
    data <- cbind(obj@frame,fitted=X %*% fixef(obj))
    ## add se bands
    data$fittedLow <- data$fitted - 1.96 * se
    data$fittedUp <- data$fitted + 1.96 * se
    if (!missing(bindAlso)) data <- cbind(data, bindAlso)
    env <- parent.frame()
    if (!missing(subset)) data <- data[eval(substitute(subset), data, env),]
    geom_ribbon(data=data, aes(ymin=fittedLow, ymax=fittedUp), alpha = alpha, ...)
}

