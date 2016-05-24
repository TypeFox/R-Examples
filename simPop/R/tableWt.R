#' Weighted cross tabulation
#' 
#' Compute contingency tables taking into account sample weights.
#' 
#' For each combination of the variables in \code{x}, the weighted number of
#' occurence is computed as the sum of the corresponding sample weights.  If
#' weights are not specified, the function \code{\link{table}} is applied.
#' 
#' @name tableWt
#' @param x a vector that can be interpreted as a factor, or a matrix or
#' \code{data.frame} whose columns can be interpreted as factors.
#' @param weights an optional numeric vector containing sample weights.
#' @param useNA a logical indicating whether to include extra \code{NA} levels
#' in the table.
#' @return The (weighted) contingency table as an object of class \code{table},
#' an array of integer values.
#' @author Andreas Alfons and Stefan Kraft
#' @seealso \code{\link{table}}, \code{\link{contingencyWt}}
#' @keywords category
#' @export
#' @examples
#' 
#' data(eusilcS)
#' tableWt(eusilcS[, c("hsize", "db040")], weights = eusilcS$rb050)
#' tableWt(eusilcS[, c("rb090", "pb220a")], weights = eusilcS$rb050, 
#'     useNA = "ifany")
#' 
tableWt <- function(x, weights = NULL, useNA = c("no", "ifany", "always")) {
    # initializations
    if(!is.data.frame(x)) x <- as.data.frame(x)
    if(is.null(weights)) return(table(x, useNA=useNA))
    else if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
    else if(length(weights) != nrow(x)) {
        stop("length of 'weights' must equal the number of rows in 'x'")
    } else if(!all(is.finite(weights))) stop("missing or infinite weights")
    useNA <- match.arg(useNA)
    if(nrow(x) > 0 && ncol(x) > 0 && useNA != "no") {
        always <- useNA == "always"
        if(ncol(x) == 1) x[, 1] <- factorNA(x[, 1], always)
        else x <- as.data.frame(lapply(x, factorNA, always))
    }
    # compute and return weighted table
    tab <- round(tapply(weights, x, sum))
    tab[is.na(tab)] <- 0
    class(tab) <- "table"
    return(tab) 
}
