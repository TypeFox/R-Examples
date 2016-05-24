#' Short \code{qselection} summary
#'
#' \code{\link{qselection}} summary
#' @param x \code{qselection} object.
#' @param \ldots Other options.
#' @return The function returns a summary table with the subsets of size
#'   \eqn{q}, their information criterion values and the chosen variables for
#'   each one. Additionally, an asterisk is shown next to the size of subset
#'   which minimizes the information criterion.
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' @seealso \code{\link{selection}}.
#' @examples
#' library(FWDselect)
#' data(diabetes)
#' x = diabetes[ ,2:11]
#' y = diabetes[ ,1]
#' obj2 = qselection(x, y, qvector = c(1:9), method = "lm", criterion = "variance", cluster = FALSE)
#' obj2
#' @export


print.qselection <- function(x = object, ...) {
  if (inherits(x, "qselection")) {
    object = x
    aux = cbind(object[[1]], object[[2]], as.character(object[[3]]))
    colnames(aux) = names(object)
    aux2 = as.data.frame(aux)
    ii <- which.min(object[[2]])
    best <- rep("", length(object[[1]]))
    best[ii] <- "*"
    res <- cbind(aux2,best)
    colnames(res) <- c(names(object), "")
    print(res)
  }else{
    stop("Argument x must be either qselection object.")
  }
}
