gauss.cor.test.numeric <- function (x, y, ...)
{
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

    if (!is.vector(x) || !is.vector(y) || length(x) != length(y))
        stop("'x' and 'y' must be vectors with the same length")

    res <- cor.test(qnorm(rank(x) / (length(x) + 1)),
                    qnorm(rank(y) / (length(y) + 1)),
                    method="pearson", ...)

    res$method <- "Gaussian rank correlation estimator"
    res$data.name <- DNAME

    res
}

setMethod("gauss.cor.test", signature(x="numeric", y="numeric"),
          gauss.cor.test.numeric)

gauss.cor.test.formula <- function(x, y, na.action, ...)
{
    if (length(x) != 2L)
        stop("formula invalid")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$y, parent.frame())))
        m$y <- as.data.frame(y)
    m[[1L]] <- as.name("model.frame")
    m$... <- NULL
    m$formula <- m$x
    m$data <- m$y
    m$x <- NULL
    m$y <- NULL
    mf <- eval(m, parent.frame())
    if (length(mf) != 2L)
        stop("formula invalid")
    DNAME <- paste(names(mf), collapse = " and ")

    ret <- gauss.cor.test(mf[[1]], mf[[2]], ...)
    ret$data.name <- DNAME
    ret
}

setMethod("gauss.cor.test", signature(x="formula", y="data.frame"),
          gauss.cor.test.formula)
