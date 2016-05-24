#' compute E^2(z) of given cosmology
#' @description Compute the \eqn{E^2(z)=(H(z)/H_0)^2}
#' @seealso \code{\link{eZ}}
#' @param z Redshift
#' @param cosmo cosmology parameter list, contains 'omegaM0', 'omegaL0', 'omegaK'
#' @return The dimentionless Hubble constant \eqn{H(z)/H_0}
#' @references Equation (14) in [H99]
#' 
#' @examples
#' eZ2(1.2,parameter.fidcosmo)
eZ2 <-
function(z, cosmo) {
  if ("w" %in% names(cosmo)) {
    return(cosmo[['omegaM0']] * (1+z)^3. +
           cosmo[['omegaL0']] * (1+z)^(1+cosmo[['w']])+
           cosmo[['omegaK']] * (1+z)^2.)
  } else {
    return(cosmo[['omegaM0']] * (1+z)^3. + cosmo[['omegaL0']] +
           cosmo[['omegaK']] * (1+z)^2.)
  }
}
