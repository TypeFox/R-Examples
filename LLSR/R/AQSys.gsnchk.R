#' @import rootSolve
require(rootSolve)
#'Merchuk's Method of calculating phase's composition in a Tieline applied to Murugesan Binodal Function
#' @rdname AQSys.gsnchk
#' @name AQSys.gsnchk
#' @title Merchuk's Method of calculating phase's composition in a Tieline applied to Murugesan Binodal Function
#' @description Merchuk et al. described a very straightforward method to calculate the concentration of each
#'component in the tieline giving only its global composition and phase's properties (such as volume and density).
#'However, other researchers relate to have achieved better fitting results using Murugesan's functions.
#'This method offers Merchul's ways o calculate Tieline's composition but using the equation proposed by Murugesan.
#' @details Using the binodal data, the global composition of a chosen tieline and its phases properties (more precisely
#'each phase density and volume). Using the data included in LLSR package the function couldn't achieve steady-state and
#'consecutively have a poor convergence tolerance. Use for your own risk.
#' @export AQSys.gsnchk
#' @export
#' @param XYdt - Standard bidimensional data.frame used in most of functions available in this package.
#' @param Xm - Component X's concentration in the tieline's global composition.
#' @param Ym - Component Y's concentration in the tieline's global composition.
#' @param Vt - Tieline's TOP phase volume.
#' @param Vb - Tieline's BOTTOM phase volume.
#' @param dyt - Tieline's TOP phase density
#' @param dyb - Tieline's BOTTOM phase density
#' @param ... Additional optional arguments. None are used at present.
#' @return sysres - The function returns the Critical Point (X,Y), Tieline Length (TLL), Tieline's Equivolume point (xVRe2o,yVRe2o),
#'and Tieline's Slope.
#' @examples
#' #
#' XYdt <- peg4kslt[,1:2]
#' #
#' Xm <- peg4kslt[2,3]
#' Ym <- peg4kslt[2,4]
#' Vt <- peg4kslt[2,5]
#' Vb <- peg4kslt[2,6]
#' dyt <- peg4kslt[2,7]
#' dyb <- peg4kslt[2,8]
#' #
#' AQSys.gsnchk(XYdt,Xm,Ym,Vt,Vb,dyt,dyb)
AQSys.gsnchk <- function(XYdt,Xm,Ym,Vt,Vb,dyt,dyb,...) {
  # Fit XYdt data to Murugesan's equation and store it in Smmry
  Smmry <- summary(mrgsn(XYdt))
  # extract regression parameters from Smmry
  A <- Smmry$coefficients[1]
  B <- Smmry$coefficients[2]
  C <- Smmry$coefficients[3]
  # Calculate alfa for a given system composition
  alfa <- Vt * dyt / (Vt * dyt + Vb * dyb)
  # the system of equations below was uses the method described by Merchuk
  # in the manuscript Merchuk, J.C., B.A. Andrews, and J.A. Asenjo. (1998),
  # Aqueous two-phase systems for protein separation: Studies on phase inversion.
  # Journal of Chromatography B: Biomedical Sciences and Applications. 711(1-2):
  #p. 285-293. to calculate tieline composition using murugesan's equation.
  # the lines below set the equation's system
  sys <- function(x) {
    F1 <- A + B * (x[2]) ^ 0.5 + C * x[2] - x[1]
    F2 <- A + B * (x[4]) ^ 0.5 + C * x[4] - x[3]
    F3 <- (Ym / alfa) - ((1 - alfa) / alfa) * x[3] - x[1]
    F4 <- (Xm / alfa) - ((1 - alfa) / alfa) * x[4] - x[2]
    #
    c(
      F1 = F1, F2 = F2, F3 = F3, F4 = F4
    )
  }
  # solve the system of equation for a given set of guess and restricting of positive
  # only results
  (sysres <- multiroot(
    f = sys, start = c(1,0,0,1),positive = TRUE
  ))
  # Calculate the tieline length and store it in sysres under the TLL alias
  sysres$TLL <-
    sqrt((sysres$root[1] - sysres$root[3]) ^ 2 + (sysres$root[2] - sysres$root[4]) ^
           2)
  # set alfa to 0.5 to calculate concentration at equivolume point
  alfaVRe2o <- 0.5
  # calculate the system composition at equivolume
  sysres$yVRe2o <-
    alfaVRe2o * (sysres$root[1] + sysres$root[3] * ((1 - alfaVRe2o) / alfaVRe2o))
  sysres$xVRe2o <-
    alfaVRe2o * (sysres$root[2] + sysres$root[4] * ((1 - alfaVRe2o) / alfaVRe2o))
  # set var name for root results (phase's composition for a given tieline)
  names(sysres$root) <- c("YT","XT","YB","XB")
  # calculate and store tieline's slope
  sysres$S <-
    (sysres$root["YT"] - sysres$root["YB"]) / (sysres$root["XT"] - sysres$root["XB"])
  # removing Slope's header to make easier its retrieve
  names(sysres$S) <- NULL
  # return all calculated parameters
  sysres
}
