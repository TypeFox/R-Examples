############ RR cloglog link function######

#' Log-Log link function with Randomized Response parameters.
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
RRlink.cloglog <- function(c,d) {
  ## link
  linkfun <- function(y) log(-log(((c+d)-y)/d))
  ## inverse link
  linkinv <- function(eta) (c+d)-d*exp(-exp(eta))
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) d*exp(eta)*exp(-exp(eta))
  valideta <- function(eta) TRUE
  link <- "RRcloglog"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta,
                 name = link),class = "link-glm")
}
