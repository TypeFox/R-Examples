############ RR cauchit link function######

#' Cauchit link function with Randomized Response parameters.
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
RRlink.cauchit <- function(c,d) {
  ## link
  linkfun <- function(y) tan(((y-c)/d - 0.5)*pi)
  ## inverse link
  linkinv <- function(eta) c + d*(0.5 + atan(eta)/pi)
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) {d/(pi*(1+eta**2))}
  valideta <- function(eta) TRUE
  link <- "RRcauchit"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta,
                 name = link),class = "link-glm")
}
