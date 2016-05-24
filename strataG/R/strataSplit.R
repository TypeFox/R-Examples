#' @title Split Strata
#' @description Return a list of \code{gtypes} for each stratum.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param strata a character vector giving a subset of strata to select. 
#'   If \code{NULL} then a list with all strata is created.
#' @param remove.sequences logical. If \code{TRUE} any sequences not referenced 
#'   in selected samples will not be in the returned object.
#' 
#' @return A named list where each element is a \code{gtypes} object 
#'   for a single stratum in \code{g}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
strataSplit <- function(g, strata = NULL, remove.sequences = FALSE) {
  if(is.null(strata)) strata <- strataNames(g)
  sapply(strata, function(st) {
    g[, , st, drop = remove.sequences]
  }, simplify = FALSE)
}