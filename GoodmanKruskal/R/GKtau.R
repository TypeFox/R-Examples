#' Compute Goodman and Kruskal tau measure of association.
#'
#' \code{GKtau} returns forward and backward Goodman and Kruskal tau measures
#' between categorical variables.
#'
#' The Goodman and Kruskal tau measure is an asymmetric association measure
#' between two categorical variables, based on the extent to which variation
#' in one variable can be explained by the other.  This function returns a
#' dataframe with both forward and backward associations.
#'
#' @param x A categorical vector (factor).
#' @param y A categorical vector (factor).
#' @param dgts Integer, number of digits for results; optional (default = 3).
#' @param includeNA Character, passed to useNA parameter for table;
#' default is "ifany"; other valid options are "no" and "always"
#' @return A one-row dataframe with the following columns:
#' \itemize{
#'   \item the names of the x and y variables,
#'   \item the numbers of distinct values Nx and Ny for each variable, and
#'   \item the forward and backward associations, tau(x,y) and tau(y,x).
#' }
#' @examples
#' x <- rep(c("a", "b", "c", "d"), each = 3)
#' y <- rep(c("a", "b", "c", "d"), times = 3)
#' z <- rep(c("a", "b", "a", "c"), each = 3)
#' GKtau(x, y)
#' GKtau(x, z)
#' GKtau(y, z)
#'
#' @author Ron Pearson
#' @export
#'
GKtau <- function(x, y, dgts = 3, includeNA = "ifany"){
  #
  #  Retrieve the names of the variables x and y
  #
  xName <- deparse(substitute(x))
  yName <- deparse(substitute(y))
  #
  #  Compute the joint empirical distribution PIij
  #
  Nij <- table(x, y, useNA = includeNA)
  PIij <- Nij/sum(Nij)
  #
  #  Compute the marginals
  #
  PIiPlus <- rowSums(PIij)
  PIPlusj <- colSums(PIij)
  #
  #  Compute marginal and conditional variations
  #
  vx <- 1 - sum(PIiPlus^2)
  vy <- 1 - sum(PIPlusj^2)
  xyTerm <- apply(PIij^2, MARGIN = 1, sum)
  vyBarx <- 1 - sum(xyTerm/PIiPlus)
  yxTerm <- apply(PIij^2, MARGIN = 2, sum)
  vxBary <- 1 - sum(yxTerm/PIPlusj)
  #
  #  Compute forward and reverse associations
  #
  tauxy <- (vy - vyBarx)/vy
  tauyx <- (vx - vxBary)/vx
  #
  #  Form summary dataframe and return
  #
  sumFrame <- data.frame(xName = xName, yName = yName,
                         Nx = nrow(Nij), Ny = ncol(Nij),
                         tauxy = round(tauxy, digits = dgts),
                         tauyx = round(tauyx, digits = dgts))
  return(sumFrame)
}
