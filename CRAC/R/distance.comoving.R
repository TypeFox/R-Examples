#' compute the comoving distance (line-of-sight) [Mpc/h]
#' 
#' @param z redshift upper limit
#' @param cosmo cosmological parameter list
#' @param z0 redshift lower limit
#' @param ... pass to integrate() to control integration accuracy.
#' @return Comoving distance from \eqn{z_0(=0)} to \eqn{z} to z [Mpc/h]
#' @references Equation (15) in [H99]
#' 
#' @examples
#' distance.comoving(0.2,parameter.fidcosmo,z0=0.3)
#' sapply(seq(0,1,0.1),function (x) distance.comoving(x,parameter.fidcosmo))
#' 
distance.comoving <-
function(z, cosmo, z0=0, ...) {
  parameter.DH * integrate(function(x){1/eZ(x,cosmo)},z0,z,...)$value
}
