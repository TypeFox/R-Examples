#' Compute the luminosity distance [Mpc/h]
#'
#' @param z redshift upper limit
#' @param cosmo cosmological parameter list
#' @param z0 redshift lower limit
#' @param ... pass to integrate() to control integration accuracy.
#' @return Luminosity distance from \eqn{z_0(=0)} to \eqn{z} [Mpc/h]
#' @references Equation (20) in [H99]
#' @seealso \code{\link{distance.angular}},\code{\link{distance.comoving}}
#' 
#' @examples
#'distance.luminosity(0.1,parameter.fidcosmo)
distance.luminosity <-
function(z, cosmo, z0=0, ...) {
  distance.angular(z, cosmo, z0=z0, ...) * (1+z)^2.
}
