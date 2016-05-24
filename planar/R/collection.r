##' Light intensity from the transmission of a bunch of plane waves at a planar interface
##'
##' Integration is performed over the solid angle defined by omega
##' @title collection_ml
##' @param xyz position matrix
##' @param wavelength wavelength
##' @param omega collection angle
##' @param psi polarisation angle
##' @param epsilon vector of permittivities
##' @param thickness thickness corresponding to each medium
##' @param maxEval passed to cubature
##' @param reqAbsError passed to cubature
##' @param tol passed to cubature
##' @param progress logical display progress bar
##' @return data.frame intensity at the x, y, z position
##' @export
##' @author Baptiste Auguie
collection_ml <- function(xyz, wavelength=632.8, omega = c(40, 50)*pi/180,
                          psi=0,
                          epsilon = c(1.5^2, epsAg(wavelength)$epsilon,
                                      1.0^2, 1.0^2),
                          thickness = c(0, 50, 10, 0),
                          maxEval = 3000, reqAbsError = 0.0,
                          tol=1e-04, progress=FALSE){

  k0 <- 2*pi/wavelength
  I <- cpp_field_collection(xyz, k0, psi, omega, epsilon, thickness,
                                   as.integer(maxEval), as.double(reqAbsError),
                                   tol, progress)

  I
}
