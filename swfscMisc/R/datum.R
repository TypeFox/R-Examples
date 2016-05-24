#' @title Datum
#' @description Return parameters specifying ellipsoid datum model.
#' 
#' @param model character, specifying which model to use for ellipsoid model. 
#'   Options are: "wgs84", "grs80", "airy", "international", "clarke", "grs67", 
#'   or partial matches thereof (case-sensitive).
#' 
#' @note Model parameters are based on distances in km.
#' 
#' @return vector of a, b, and f parameters.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
datum <- function(model = c("wgs84", "grs80", "airy", "international", "clarke", "grs67")) {
  switch(match.arg(model),
    wgs84 = c(a = 6378137, b = 6356752.3142, f = 1 / 298.257223563),
    grs80 = c(a = 6378137, b = 6356752.3141, f = 1 / 298.257222101),
    airy = c(a = 6377563.396, b = 6356256.909, f = 1 / 299.3249646),
    international = c(a = 6378888, b = 6356911.946, f = 1 / 297),
    clarke = c(a = 6378249.145, b = 6356514.86955, f = 1 / 293.465),
    grs67 = c(a = 6378160, b = 6356774.719, f = 1 / 298.25)
  )
}
