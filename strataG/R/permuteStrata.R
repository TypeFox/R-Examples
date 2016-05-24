#' @title Permute strata
#' @description Permute the strata slot within a \linkS4class{gtypes} object.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return a \linkS4class{gtypes} object with the strata randomly permuted.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' summary(msats.g)
#' 
#' ran.msats <- permuteStrata(msats.g)
#' summary(ran.msats)
#'
#' @export
#' 
permuteStrata <- function(g) {
  st <- strata(g)
  no.nas <- st[!is.na(st)]
  no.na.sample <- sample(as.character(no.nas))
  names(no.na.sample) <- names(no.nas)
  strata(g) <- factor(no.na.sample[names(st)])
  g
}