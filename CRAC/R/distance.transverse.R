#' compute the comoving distance (transverse) [Mpc/h]
#' 
#' @param z redshift upper limit
#' @param cosmo cosmological parameter list
#' @param z0 redshift lower limit
#' @param ... pass to integrate() to control integration accuracy.
#' @return Comoving distance from \eqn{z_0(=0)} to \eqn{z} [Mpc/h]
#' @seealso \code{\link{distance.comoving}}
#' @references Equation (16) in [H99]
#' 
#' @examples
#' distance.transverse(0.1,parameter.fidcosmo)
distance.transverse <-
function(z, cosmo, z0=0, ...) {
  dc <- distance.comoving(z, cosmo, z0=z0)
  if (cosmo[['omegaK']] > 0){
    parameter.DH * sinh(sqrt(cosmo[['omegaK']])*dc/parameter.DH) /
        sqrt(cosmo[['omegaK']])
  } else if (cosmo[['omegaK']] < 0){
    parameter.DH * sin(sqrt(-cosmo[['omegaK']])*dc/parameter.DH) /
        sqrt(-cosmo[['omegaK']]) 
  } else {
    dc
  }
}
