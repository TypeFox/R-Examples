#' compute the angular diameter distance [Mpc/h]
#' 
#' @description ONLY FOR flat universe, \eqn{\Omega_\mathbf{k}=0}.
#' @param z redshift upper limit
#' @param cosmo cosmological parameter list
#' @param z0 redshift lower limit
#' @param ... pass to integrate() to control integration accuracy.
#' @return Angular diameter distance from \eqn{z_0(=0)} to \eqn{z} [Mpc/h]
#' @references Equation (18) in [H99]  
#' @seealso \code{\link{distance.comoving}}
#' 
#' @examples
#' distance.angular(0.1,parameter.fidcosmo)
#' distance.angular(0.3,list(omegaM0=0.272,omegaL0=0.728,omegaK=0.0,h=0.704))
#' 
distance.angular <-
function(z, cosmo, z0=0, ...) {
  dc1 <- distance.comoving(z, cosmo, ...)
  if (z0 == 0) {
    dc1/(1+z)
  } else {
    dc2 <- distance.comoving(z0,cosmo, ...)
    (dc1-dc2)/(1+z0) #' !! ONLLY for omegaK=0
  }
}
