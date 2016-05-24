#' calculate the critical density at redshift z
#'
#' @param z redshift
#' @param cosmo cosmological parameter list
#'
#' return the critical densith at given redshift in \eqn{h_0^2 \mathrm{kg/m}^3]}
#' 
#' @examples
#'  # get the critial density at z=0
#'  rhoc(0,parameter.fidcosmo) 
#' 
rhoc <-
function(z,cosmo) {
  rhoc0 <- 3*(1e5/3.0856776e22)^2/8./pi/parameter.constant[['G']]
  return (rhoc0*eZ2(z,cosmo))
}
