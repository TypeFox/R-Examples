#'  @title Get units of a scanonevar object
#'
#'  @author Robert Corty \email{rcorty@@gmail.com}
#'
#'  @param object the scanonevar object whose units are interrogated
#'
#'  @description Utility function to get the units of a scanonevar, which is either 'LODs'
#'    (logarithm of the odds between null and alternative model) or 'emp.ps' (empirically
#'    determined p.value of that LOD)
#'
#'  @return units of scanonevar object
#'
#'  @details none
units <- function(object) {

  return(attr(object, 'units'))
}
