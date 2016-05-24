#' Compute the norm of a vector
#' @description Computes the norm of a vector with size adjustment.
#' @param x A numeric vector.
#' @param norm The type of norm. See details.
#' @details This function calculates the norm of a vector. It is slightly different from the
#' ordinary norm function in that \code{L1}-norm and \code{L2}-norm are adjusted by the length
#' of vector, \eqn{n}. To be specific, define a vector \eqn{x = (x_1,x_2,...,x_n)}.
#' The \code{L1}-norm is defined as
#' \deqn{||x||_1 = (|x_1| + |x_2| + ... + |x_n|)/n.}
#' The \code{L2}-norm is defined as
#' \deqn{||x||_2 = (\sqrt{(x_1)^2 + (x_2)^2 + ... + (x_n)^2)/n.}}
#' The \code{Linf}-norm is defined as
#' \deqn{||x||_{\infty} = \max(|x_1|,|x_2|,...,|x_n|).}
#' Note that a matrix \code{X} will be coerced to a vector.
#' @return A non-negative number.
#'
#' @author Debin Qiu, Jeongyoun Ahn
#' @examples x <- 1:10
#' # "L1" norm, same as norm(as.matrix(x), "1")*length(x)
#' norm_vec(x)
#'
#' # "L2" norm
#' norm_vec(x, "L2")
#'
#' # "Linf" norm, same as norm(as.matrix(x),"I")
#' norm_vec(x, "Linf")
#'
#' @export
#'
norm_vec <- function(x,norm = c("L1","L2","Linf")) {
  norm <- match.arg(norm)
  switch(norm, L1 = sum(abs(x))/length(x),
         L2 = sum(x^2)/length(x), Linf = max(abs(x)))
}
