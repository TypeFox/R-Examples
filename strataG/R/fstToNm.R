#' @title Nm from Fst
#' @description Calculate Nm (number of migrants per generation) for a 
#'   given value of Fst.
#'
#' @param fst estimate of Fst between populations.
#' @param ploidy ploidy of locus Fst is from.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#' 
fstToNm <- function(fst, ploidy) {
  ((1 / fst) - 1) / (ploidy * 2)
}