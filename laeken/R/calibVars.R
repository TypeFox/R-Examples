# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

#' Construct a matrix of binary variables for calibration
#'
#' Construct a matrix of binary variables for calibration of sample weights
#' according to known marginal population totals.
#'
#' @name calibVars
#' @param x a vector that can be interpreted as factor, or a matrix or
#' \code{data.frame} consisting of such variables.
#'
#' @return A matrix of binary variables that indicate membership to the
#' corresponding factor levels.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link{calibWeights}}
#'
#' @keywords survey
#'
#' @examples
#' data(eusilc)
#' # default method
#' aux <- calibVars(eusilc$rb090)
#' head(aux)
#' # data.frame method
#' aux <- calibVars(eusilc[, c("db040", "rb090")])
#' head(aux)
#'
#' @export

calibVars <- function(x) UseMethod("calibVars")

#' @export
calibVars.default <- function(x) {
    if(length(x) == 0) matrix(integer(), 0, 0)
    x <- as.factor(x)
    res <- sapply(levels(x), function(l) as.integer(x == l))
    rownames(res) <- names(x)  # set rownames from original vector
    res
}

#' @export
calibVars.matrix <- function(x) calibVars(as.data.frame(x))

#' @export
calibVars.data.frame <- function(x) {
    res <- lapply(x, calibVars)  # list of matrices for each variable
    res <- mapply(function(x, nam) {
            colnames(x) <- paste(nam, colnames(x), sep=".")
            x
        }, res, names(x), SIMPLIFY=FALSE)
    res <- do.call("cbind", res)  # combine matrices
    rownames(res) <- row.names(x)  # set rownames from original data.frame
    res
}
