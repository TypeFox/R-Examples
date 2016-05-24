

#' @rdname AQSys.tielines
#' @title Merchuk's Method - Tieline's Composition Calculation
#' @description Merchuk et al. described a very straightforward method to calculate the concentration of each component in the
#' tieline giving only its global composition and phase's properties (such as volume and density).
#' @details Using the binodal data, the global composition of a chosen tieline and its phases properties (more precisely each
#' phase density and volume)
#' @method AQSys tielines
#' @export AQSys.tielines
#' @export
#' @param XYdt - Binodal Experimental data that will be used in the nonlinear fit
#' @param Xm - Component X's concentration in the tieline's global composition.
#' @param Ym - Component Y's concentration in the tieline's global composition.
#' @param Vt - Tieline's TOP phase volume.
#' @param Vb - Tieline's BOTTOM phase volume.
#' @param dyt - Tieline's TOP phase density
#' @param dyb - Tieline's BOTTOM phase density
#' @param ... Additional optional arguments. None are used at present.
#' @return The function returns the Critical Point (X,Y), Tieline Length (TLL), Tieline's Equivolume point (xVRe2o,yVRe2o),
#' and Tieline's Slope.
#' @examples
#' \dontrun{
#' AQSys.tielines(XYdt,Xm,Ym,Vt,Vb,dyt,dyb)
#' }
AQSys.tielines <- function(XYdt,Xm,Ym,Vt,Vb,dyt,dyb,...) {
  # Fit XYdt data to Merchuk's equation and store it in Smmry
  Smmry <- summary(mrchk(XYdt))
  # extract regression parameters from Smmry
  P1 <- Smmry$coefficients[1]
  P2 <- Smmry$coefficients[2]
  P3 <- Smmry$coefficients[3]
  # Calculate alfa for a given system composition
  alfa <- Vt * dyt / (Vt * dyt + Vb * dyb)
  # the system of equations below was uses the method described by Merchuk
  # in the manuscript Merchuk, J.C., B.A. Andrews, and J.A. Asenjo. (1998),
  # Aqueous two-phase systems for protein separation: Studies on phase inversion.
  # Journal of Chromatography B: Biomedical Sciences and Applications. 711(1-2):
  #p. 285-293. to calculate tieline composition using merchuk's equation.
  # the lines below set the equation's system
  sys <- function(x) {
    F1 <- P1 * exp(P2 * x[2] ^ 0.5 - P3 * x[2] ^ 3) - x[1]
    F2 <- P1 * exp(P2 * x[4] ^ 0.5 - P3 * x[4] ^ 3) - x[3]
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
#' @rdname AQSys.crpt
#' @title Merchuk's Method - Critical Point Calculation
#' @description Merchuk et al. described a very straightforward method to calculate the critical composition of a given Binodal
#' curve and its Tielines.
#' @details Using the binodal data, tieline's Slopes (S), the composition of bottom-rich component in the bottom phase (XB)
#' and an equation which stablish a relatioship between them, this function returns the critical composition of the binodal
#' under study. When used within a iterative function, mrchk.tielines() can be used to obtain TLL and S and therefore
#' calculate the critical composition.
#' @method AQSys crpt
#' @export AQSys.crpt
#' @export
#' @param tldata - A data.frame with two columns containing a set of Tieline's Slopes (S)
#' and its bottom-rich component composition in the bottom phase (XB).
#' @param XYdt - Binodal Experimental data that will be used in the nonlinear fit
#' @param ... Additional optional arguments. None are used at present.
#' @return (XCrit,YCrit) - The function returns Tieline's Critical Point Composition
#' @examples
#' \dontrun{
#' AQSys.crpt(XYdt,tldata)
#' }
AQSys.crpt <- function(XYdt, tldata,...) {
  # Fit XYdt data to Merchuk's equation and store it in Smmry
  Smmry <- summary(mrchk(XYdt))
  # store tieline data into a dataframe variable. It might be a better approach check if
  # user stored it in a dataframe and if not trigger an error.
  tldata <- as.data.frame(tldata)
  # calculate the regression coefficient of a polynomial equation of order 2
  resTLLfit <- lm(S ~ poly(XB, 2, raw = TRUE), data = tldata)
  # solve a system of equations using Merchuk's equation and the
  # polynomial from tieline regression
  fitC <- function(x) {
    F1 = -x[1] + Smmry$coefficients[1] * exp(Smmry$coefficients[2] * x[2] ^
                                               0.5 - Smmry$coefficients[3] * x[2] ^ 3) * ((Smmry$coefficients[2] / (2 *
                                                                                                                      x[2] ^ 0.5)) - 3 * Smmry$coefficients[3] * x[2] ^ 2)
    F2 = -x[1] + resTLLfit$coefficients[1] + resTLLfit$coefficients[2] *
      x[2] + resTLLfit$coefficients[3] * x[2] ^ 2
    c(F1 = F1,F2 = F2)
  }
  # solve the system of equation for a given set of guess
  (sres <- multiroot(f = fitC, start = c(.1,.1)))
  # store calculated heavy-component concentration at global composition
  XCrit <- sres$root[2]
  # store calculated light-component concentration at global composition
  YCrit <-
    Smmry$coefficients[1] * exp(Smmry$coefficients[2] * (XCrit ^ (0.5)) - Smmry$coefficients[3] *
                                  (XCrit ^ 3))
  # include critical concentration in the output of the system of equations result
  sres$XCrit <- XCrit
  sres$YCrit <- YCrit
  # return all calculated parameters
  sres
}
