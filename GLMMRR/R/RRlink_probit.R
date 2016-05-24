############ RR probit link function######

#' Probit link function with Randomized Response parameters.
#'
#' Reference: Fox, J-P, Klotzke, K. and Veen, D. (2016). \emph{Generalized Linear Mixed Models for Randomized
#' Responses.} Manuscript submitted for publication.
#'
#' @param c
#' a numeric vector containing the parameter c.
#' @param d
#' a numeric vector containing the parameter d.
#' @return
#' RR link function.
#'
#' @export
RRlink.probit <- function(c,d) {
  ## link
  linkfun <- function(y) qnorm((y-c)/d)
  ## inverse link
  linkinv <- function(eta) c + d*pnorm(eta)
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) {d*dnorm(eta)}
  valideta <- function(eta) TRUE
  link <- "RRprobit"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta,
                 name = link),class = "link-glm")
}
