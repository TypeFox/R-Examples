##' @title Compare proxies across two data sets
##'
##' @description.. content for \description{} (no empty lines) ..
##'
##' @details To Do
##'
##' @param x data frame; training set samples to compare against
##' @param ... other arguments passed to methods.
##'
##' @return To Do
##'
##' @author Gavin L. Simpson
##'
##' @rdname compare
##'
 `compare` <- function(x, ...) {
    UseMethod("compare")
}

##' @param y data frame; passive or core samples
##' @param env numeric vector of environmental or contraint data for residual length ordination
##' @param by character; compare data sets by sites or species (proxies).
##' @param ordination character; which constrained ordination method to use
##' @param method character; which dissimilarity method to use. See \code{distance}.
##' @param transform character; should a transformation be applied to the data. Ignored.
##' @param n2limit integer; the criterion for indicating species with potentially
##' poorly estimated optima. The default value of \code{5L} is one suggested by
##' R. Telford.
##'
##' @rdname compare
##'
`compare.default` <- function(x, y, env,
                              by = c("sites", "species"),
                              ordination = "rda",
                              method = "chord",
                              transform = NULL,
                              n2limit = 5L,
                              ...) {
    joint <- join(x, y, split = TRUE)
    x2 <- joint$x
    y2 <- joint$y

    n2.x <- n2(x2, which = "species")

    by <- match.arg(by)

    if (identical(by, "species")) {
        n2.y <- n2(y2, which = "species")
        cny <- colnames(y)
        cnx <- colnames(x)
        out <- data.frame(inTrain  = cny %in% cnx,
                          n2       = n2.y[cny],
                          n2Train  = n2.x[cny],
                          max      = apply(y2, 2, max)[cny],
                          maxTrain = apply(x2, 2, max)[cny])
    } else {
        out <-
            data.frame(sumMissing  = rowSums(y2[, n2.x == 0L]),
                       sumPoorOpt  = rowSums(y2[, n2.x <= n2limit]),
                       closestSamp = minDC(ana.fit <- analog(x2, y2))$minDC)
        if (!missing(env)) {
            out <- cbind(out,
                         residLen = residLen(x, env, y,
                                             method = ordination)[["passive"]])
        }
    }
    out
}
