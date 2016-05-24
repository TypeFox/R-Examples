#' Visualization of \code{qselection} object
#'
#' This function plots the cross-validation information criterion for several
#' subsets of size \eqn{q} chosen by the user.
#' @param x \code{qselection} object.
#' @param y NULL
#' @param ylab NULL
#' @param \ldots Other options.
#' @return Simply returns a plot.
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' @seealso \code{\link{selection}}.
#' @examples
#' library(FWDselect)
#' data(diabetes)
#' x = diabetes[ ,2:11]
#' y = diabetes[ ,1]
#' obj2 = qselection(x, y, qvector = c(1:9), method = "lm", criterion = "variance", cluster = FALSE)
#' plot(obj2)
#' @importFrom graphics plot
#' @importFrom graphics axis
#' @export

plot.qselection <- function(x = object, y = NULL,
    ylab = NULL, ...) {
    object = x
    x = object$q
    y = object[[2]]
    if (is.null(ylab)) {
        ylab_aux = names(object)[2]
        if (ylab_aux == "deviance")
            ylab = "Deviance"
        if (ylab_aux == "variance")
            ylab = "Residual variance"
        if (ylab_aux == "R2")
            ylab = expression(R^2)
    }
    plot(x, y, type = "o", bg = "black", pch = 21,
        xlab = "Subset size q", ylab = ylab, xaxt = "n",
        ...)
    axis(1, x, ...)
}
