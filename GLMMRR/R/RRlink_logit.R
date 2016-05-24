############ RR logit link function######

#' Logit link function with Randomized Response parameters.
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
RRlink.logit <- function(c,d) {
  ## link
  linkfun <- function(y) log((y-c)/(c+d-y))
  ## inverse link
  linkinv <- function(eta) c + d*(exp(eta)/(1+exp(eta)))
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) {d*exp(eta)/((1+exp(eta))**2)}
  valideta <- function(eta) TRUE
  link <- "RRlogit"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta,
                 name = link),class = "link-glm")
}
