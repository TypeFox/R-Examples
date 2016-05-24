#' Weighted contingency coefficients
#'
#' Compute (weighted) pairwise contingency coefficients.
#'
#' The function \code{\link{tableWt}} is used for the computation of the
#' corresponding pairwise contingency tables. The following methods are implemented:
#' \itemize{
#' \item \code{contingencyWt.default(x, y, weights = NULL, ...)}
#' \item \code{contingencyWt.matrix(x, weights = NULL, ...)}
#' \item \code{contingencyWt.data.frame(x, weights = NULL, ...)}
#' }
#' Additional parameters are:
#' \itemize{
#' \item y: a vector that can be interpreted as factor (for the default method)
#' \item weights: an optional numeric vector containing sample weights
#' }
#'
#' @name contingencyWt
#' @param x for the default method, a vector that can be interpreted as factor.
#' For the matrix and \code{data.frame} methods, the columns should be
#' interpretable as factors.
#' @param \dots for the generic function, arguments to be passed down to the
#' methods, otherwise ignored.
#' @return For the default method, the (weighted) contingency coefficient of
#' \code{x} and \code{y}.
#'
#' For the matrix and \code{data.frame} method, a matrix of (weighted) pairwise
#' contingency coefficients for all combinations of columns.  Elements below
#' the diagonal are \code{NA}.
#' @author Andreas Alfons and Stefan Kraft
#' @export
#' @keywords methods
#' @seealso \code{\link{tableWt}}
#' @references Kendall, M.G. and Stuart, A. (1967) \emph{The Advanced Theory of
#' Statistics, Volume 2: Inference and Relationship}. Charles Griffin & Co Ltd,
#' London, 2nd edition.
#' @keywords category
#' @examples
#'
#' data(eusilcS)
#'
#' ## default method
#' contingencyWt(eusilcS$pl030, eusilcS$pb220a, weights = eusilcS$rb050)
#'
#' ## data.frame method
#' basic <- c("age", "rb090", "hsize", "pl030", "pb220a")
#' contingencyWt(eusilcS[, basic], weights = eusilcS$rb050)
contingencyWt <- function(x, ...) UseMethod("contingencyWt")

#' @export
contingencyWt.default <- function(x, y, weights = NULL, ...) {
    tab <- tableWt(data.frame(x, y), weights)
    tab <- tab[rowSums(tab) > 0, colSums(tab) > 0]
    chisq <- as.numeric(chisq.test(tab)$statistic)
    return(sqrt(chisq / (sum(tab) + chisq)))
}

#' @export
contingencyWt.matrix <- function(x, weights = NULL, ...) {
    contingencyWt(as.data.frame(x), weights=weights, ...)
}

#' @export
contingencyWt.data.frame <- function(x, weights = NULL, ...) {
    # computes *pairwise* contingency coefficients
    p <- ncol(x)
    res <- matrix(NA, ncol=p-1, nrow=p-1)
    for(i in 1:(p-1)) {
        for(j in (i+1):p) {
            res[i, j-1] <- contingencyWt(x[, i], x[, j], weights)
        }
    }
    nam <- names(x)
    dimnames(res) <- list(nam[1:(p-1)], nam[2:p])
    res
}
