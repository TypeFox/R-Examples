#' Hypothesis Tests
#'
#' This function uses \code{\link{htest.short}}, to extract statistic and p-value from \code{htest}-classed object. Main advantage of using \code{htest} is that it's vectorised, and can accept multiple methods.
#'
#' Default parameters are read from \code{options}:
#'
#' \itemize{
#'     \item 'rapport.use.labels'.
#' }
#' @param x arguments to be passed to function specified in \code{test}
#' @param ... additional arguments for function specified in \code{test}
#' @param use.labels a logical value indicating whether variable labels should be placed in row names. If set to \code{FALSE}, output of \code{deparse(substitute(x))} will be used.
#' @param use.method.names use the string provided in \code{method} attribute of \code{htest} object
#' @param colnames a character string containing column names
#' @return a \code{data.frame} with applied tests in rows, and their results (statistic and p-value) in columns
#' @examples \dontrun{
#' library(nortest)
#' htest(rnorm(100), shapiro.test)
#' htest(rnorm(100), lillie.test, ad.test, shapiro.test)
#' htest(mtcars, lillie.test)
#' htest(mtcars, lillie.test, ad.test, shapiro.test)
#' }
#' @export
htest <- function(x, ..., use.labels = getOption('rapport.use.labels'), use.method.names = TRUE, colnames = c('Method', 'Statistic', 'p-value')){

    test <- list(...)
    test.len <- length(test)
    test.name <- sapply(substitute(list(...)), deparse)[-1]

    if (is.atomic(x) || is.formula(x)){
        if (test.len == 1){
            res <- htest.short(each(test[[1]])(x))
            method.name <- attr(res, 'method')
        } else {
            res <- data.frame(lapply(each(test)(x), htest.short))
            method.name <- sapply(res, attr, which = 'method')
        }
        res <- data.frame(t(res))
        if (is.formula(x))
            x.nms <- deparse(substitute(x))
        else
            x.nms <- if (use.labels) label(x) else name(x)
        x.len <- 1
    } else {
        if (test.len == 1){
            res <- data.frame(lapply(lapply(x, test[[1]]), htest.short))
            method.name <- sapply(res, attr, which = 'method')
        } else {
            res <- lapply(x, function(y) lapply(each(test)(y), htest.short))
            method.name <- sapply(res[[1]], attr, which = "method")
        }
        res <- t(data.frame(res))
        x.nms <- if (use.labels) label(x) else name(x)
        x.len <- length(x)
    }

    if (use.method.names)
        test.name <- method.name

    if (nrow(res) == length(test.name))
        rn <- test.name
    else
        rn <- sprintf("%s (%s)", rep(test.name, x.len), rep(x.nms, each = test.len))

    names(rn) <- NULL
    rownames(res) <- NULL
    res <- cbind(rn, res)
    colnames(res) <- colnames

    return(res)
}


#' Extract Values from \code{htest} Objects
#'
#' Extract value of statistic and its p-value from \code{htest} object.
#' @param x \code{htest}-class object
#' @return named numeric vector with the value of statistic and its p-value
#' @examples \dontrun{
#' htest.short(shapiro.test(rnorm(100)))
#' }
#' @export
htest.short <- function(x){
    stopifnot(inherits(x, 'htest'))
    structure(c(x$statistic, p = x$p.value), method = x$method)
}
