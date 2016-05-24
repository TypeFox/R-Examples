#' Descriptive Statistics
#'
#' This function operates only on vectors or their subsets, by calculating a descriptive statistic specified in \code{fn} argument.
#' @param x a numeric variable to be summarised
#' @param subset an expression that evaluates to logical vector (defaults to \code{NULL}, in which case the function specified in \code{fun} is applied on a vector)
#' @param fn a function or a function name to be applied on a variable or it's subset
#' @param na.rm a logical value indicating whether \code{NA}'s should be removed (defaults to \code{TRUE})
#' @param ... additional arguments for function specified in \code{fn}
#' @return a numeric
univar <- function(x, subset = NULL, fn, na.rm = TRUE, ...) {

    if (base::missing(x))
        stop('Variable not specified.')

    ## subset the data
    if (!is.null(subset)){
        x.subset <- subset.default(x, subset)
        ## check if split was successful
        if (is.null(x.subset))
            warning('Data subset error, using whole dataset.')
        else
            x <- x.subset
    }

    res <- do.call(fn, list(x, na.rm = na.rm, ...))

    if (length(res) > 1)
        warning("Resulting statistic has more than one value. This shouldn't have happened!")

    return(res)
}


#' Number of Cases
#'
#' Returns the number of valid (non-\code{NA}) values in a variable. This is a wrapper around \code{\link{univar}} function with \code{\link{length}} function passed in \code{fn} argument, but with missing values previously removed. However, it's not possible to cancel \code{NA} omission with this function (doing so will yield error).
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with number of valid (non-NA) vector elements
#' @export
#' @aliases n rp.n
n <- function(...)
    univar(..., fn = function(...) length(..1))
#' @export
rp.n <- n


#' Number of Valid Cases
#'
#' Returns the number of valid (non-\code{NA}) values in a variable. This is a wrapper around \code{\link{univar}} function with \code{\link{length}} function passed in \code{fn} argument, but with missing values previously removed. However, it's not possible to cancel \code{NA} omission with this function (doing so will yield error).
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with number of valid (non-NA) vector elements
#' @export
#' @aliases nvalid rp.valid
nvalid <- function(...)
    univar(..., fn = function(...) length(na.omit(..1)))
#' @export
rp.valid <- nvalid


#' Number of Missing Cases
#'
#' Returns a number of missing (\code{NA}) values in a variable. This is a wrapper around \code{\link{univar}} function with anonymous function passed to count number of \code{NA} elements in a variable.
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with number of missing vector elements
#' @export
#' @aliases nmissing rp.missing
nmissing <- function(...)
    univar(..., fn = function(...) base::sum(is.na(..1)))
#' @export
rp.missing <- nmissing


#' Percent
#'
#' Calculates percentage of cases for provided variable and criteria specified in \code{subset} argument. Function accepts numeric, factor and logical variables for \code{x} parameter. If numeric and/or factor is provided, subsetting can be achieved via \code{subset} argument. Depending on value of \code{na.rm} argument, either valid (\code{na.rm = TRUE}) or all cases (\code{na.rm = FALSE}) are taken into account. By passing logical variable to \code{x}, a sum of (\code{TRUE}) elements is calculated instead, and valid percents are used (\code{NA} are excluded).
#' @param x a numeric variable to be summarised
#' @param subset an expression that evaluates to logical vector (defaults to \code{NULL})
#' @param na.rm should missing values be
#' @param pct print percent string too?
#' @param ... additional arguments for \code{\link{pct}} function
#' @return a numeric or string depending on the value of \code{pct}
#' @examples \dontrun{
#' set.seed(0)
#' x <- sample(5, 100, replace = TRUE)
#' percent(x > 2)
#' }
#' @export
#' @aliases percent rp.percent
percent <- function(x, subset = NULL, na.rm = TRUE, pct = FALSE, ...){
    if (is.logical(x)){
        res <- base::sum(x, na.rm = na.rm) / ifelse(na.rm, nvalid(x), length(x)) * 100
    } else {
        if (na.rm)
            res <- nvalid(x, subset) / nvalid(x) * 100
        else
            res <- nvalid(x, subset) / length(x) * 100
    }
    return (ifelse(pct, pct(res, ...), res))
}
#' @export
rp.percent <- percent


#' Minimum
#'
#' Returns the minimum of all values in a vector by passing \{code{\link{min}} as \code{fn} argument to \code{\link{univar}} function.
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with minimum value
#' @export
#' @aliases min rp.min
min <- function(...)
    univar(..., fn = base::min)
#' @export
rp.min <- min


#' Maximum
#'
#' Returns the maximum of all values in a vector by passing \{code{\link{max}} as \code{fn} argument to \code{\link{univar}} function.
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with maximum value
#' @export
#' @aliases max rp.max
max <- function(...)
    univar(..., fn = base::max)
#' @export
rp.max <- max


#' Range
#'
#' Calculates difference between the largest and the smallest value in a vector. See \code{\link{univar}} for details.
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with calculated range
#' @export
#' @aliases range rp.range
range <- function(...)
    univar(..., fn = function(...) diff(base::range(..1, ...)))
#' @export
rp.range <- range


#' Sum
#'
#' Returns the sum of variable's elements, by passing \code{\link{sum}} as \code{fn} argument to \code{\link{univar}} function.
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with sum of vector elements
#' @export
#' @aliases sum rp.sum
sum <- function(...)
    univar(..., fn = base::sum)
#' @export
rp.sum <- sum


#' Mean
#'
#' Calculates mean of given variable by passing \code{\link{sum}} as \code{fn} argument to \code{\link{univar}} function.
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with variable's mean
#' @export
#' @aliases mean rp.mean
mean <- function(...)
    univar(..., fn = mean.default)
#' @export
rp.mean <- mean


#' Standard Error of Mean
#'
#' Calculates standard error of mean for given variable. See \code{\link{univar}} for details.
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with standard error of mean
#' @export
#' @aliases se.mean rp.se.mean
se.mean<- function(...)
    sqrt(var(..., na.rm = TRUE) / nvalid(...))
#' @export
rp.se.mean <- se.mean


#' Standard Deviation
#'
#' Calculates standard deviation of given variable. See \code{\link{univar}} for details.
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with variable's standard deviation
#' @export
#' @aliases sd rp.sd
sd <- function(...)
    univar(..., fn = stats::sd)
#' @export
rp.sd <- sd


#' Variance
#'
#' Calculates variance of given variable. See \code{\link{univar}} for details.
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with variable's variance
#' @export
#' @aliases var rp.var
var <- function(...)
    univar(..., fn = stats::var)
#' @export
rp.var <- var


#' Median
#'
#' Calculates median of given variable. See \code{\link{univar}} for details.
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with variable's median
#' @export
#' @aliases median rp.median
median <- function(...)
    univar(..., fn = stats::median.default)
#' @export
rp.median <- median


#' Interquartile Range
#'
#' Calculates interquartile range of given variable. See \code{\link{univar}} for details.
#' @param ... parameters to be passed to \code{univar} function
#' @return a numeric value with variable's interquartile range
#' @export
#' @aliases iqr IQR rp.iqr
iqr <- function(...)
    univar(..., fn = stats::IQR)
#' @export
IQR <- iqr
#' @export
rp.iqr <- iqr
