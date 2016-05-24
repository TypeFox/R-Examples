#' @name pin
#' @title Title.
#' @description Description.
#' @param object Output of either asreml or the function AI_algorithm
#' @param transform The transformation used in the delta-method
#' @return blabla
#' @author VSNI?
#' @details blabla.
#' @examples blabla.
#' @export
#' @references \\url{http://homepages.ed.ac.uk/imsw/asreml/useofpin.pdf}
#'             Most of it is borrowed from Venables and Ripley, S Programming, p170
#
pin <- function (object, transform) {
  pframe <- as.list(object$gammas)
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)),pframe)
  X <- as.vector(attr(tvalue, "gradient"))
  X[object$gammas.type == 1] <- 0
  tname <- if (length(transform) == 3)
  transform[[2]]
  else ""
  n <- length(pframe)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  k <- 1 + (i > j)
  Vmat <- object$ai
  se <- sqrt(sum(Vmat * X[i] * X[j] * k))
  data.frame(row.names = tname, Estimate = tvalue, SE = se)
}