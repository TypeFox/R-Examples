#' validate the cosmologial parameters are enough
#'
#' @description validate the cosmological parameter list is complete for other usage. 
#'  If missing data detected, it will be filled with fiducial value.
#'  Also the warning message is printed in red for *nix enviroment.
#' @seealso \code{\link{parameter.fidcosmo}}
#' @param cosmo cosmological parameter list
#' @return the cosmo list will be updated if important variables are missing
#' 
#' @examples
#' # there are two typos in the cosmology parameter list
#' a <- list(omegaM=0.272,omegaL0=0.728,omegaK=0.0,h0=0.704)
#' valid.cosmo(a)
#' print(a)
valid.cosmo <-
function(cosmo) {
  if (!('omegaM0' %in% names(cosmo))) {
    cat("\033[31m Missing omegaM0 \033[0m\n")
  }
  if (!('omegaL0' %in% names(cosmo))) {
    cat("\033[31m Missing omegaL0 \033[0m\n")
  }
  if (!('omegaK' %in% names(cosmo))) {
    cat("\033[31m No k found, computed \n \033[0m")
    temp = cosmo
    temp['omegaK'] <- 1-cosmo[['omegaM0']]-cosmo[['omegaL0']]
    eval.parent(substitute(cosmo <- temp))
  }
  if (!('h' %in% names(cosmo))) {
    cat("\033[31m No h found, assigned 0.7 \n \033[0m")
    temp = cosmo
    temp['h'] <- parameter.fidcosmo[['h']]
    eval.parent(substitute(cosmo <- temp))
  }
}
