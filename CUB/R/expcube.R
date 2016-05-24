#' @title Expectation of CUBE models
#' @description Compute the expectation of a CUBE model without covariates.
#' @aliases expcube
#' @usage expcube(m, pai, csi, phi)
#' @param m Number of ordinal categories
#' @param pai Uncertainty parameter
#' @param csi Feeling parameter
#' @param phi Overdispersion parameter
#' @export expcube
#' @seealso \code{\link{varcube}}, \code{\link{varcub00}}, \code{\link{expcub00}} 
#' @keywords distribution
#' @references Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data,
#'  \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 \cr
#' Iannario, M. (2015). Detecting latent components in ordinal data with overdispersion by means 
#' of a mixture distribution, \emph{Quality & Quantity}, \bold{49}, 977--987
#' @examples
#' m<-10
#' pai<-0.1
#' csi<-0.7
#' phi<-0.2
#' meancube<-expcube(m,pai,csi,phi)


expcube <-
function(m,pai,csi,phi){expcub00(m,pai,csi)}
