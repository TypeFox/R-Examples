#' Glacier locations
#'
#' Downloaded by Dianne Cook.
#'
#' @format
#' \describe{
#' \item{id}{A 12-character unique glacier identifier}
#' \item{name}{30-character name of the glacier.}
#' \item{lat,long}{Location of the glacier: "The point on the glacier whose
#'   coordinates are given should be in the upper part of the ablation area,
#'   in the main stream and sufficiently high so as not to be lost if the
#'   glacier retreats"}
#' \item{area}{The total area of the glacier in a horizontal projection in
#'   square kilometers, up to 6 digits.}
#' \item{country}{2-character abbreviation for the name of the country or
#'   territory in which the glacier is located. These codes are ISO3166 country
#'   codes}
#' }
#'
#' @source \url{http://nsidc.org/data/glacier_inventory/}
"glaciers"
