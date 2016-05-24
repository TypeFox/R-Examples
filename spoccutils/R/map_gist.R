#' Make an interactive map to view in the browser as a Github gist
#'
#' @export
#' @importFrom gistr gist_create
#'
#' @param data A data.frame, with any number of columns, but with at least the
#'    following: name (the taxonomic name), latitude (in dec. deg.), longitude
#'    (in dec. deg.)
#' @param description Description for the Github gist, or leave to default (=no description)
#' @param public (logical) Whether gist is public (default: TRUE)
#' @param browse If TRUE (default) the map opens in your default browser.
#' @param ... Further arguments passed on to \code{\link{style_geojson}}
#'
#' @details See \code{\link[gistr]{gist_auth}} for help on authentication
#'
#' @examples \dontrun{
#' library("spocc")
#' spp <- c('Danaus plexippus','Accipiter striatus','Pinus contorta')
#' dat <- occ(spp, from=c('gbif','ecoengine'), limit=30, gbifopts=list(hasCoordinate=TRUE))
#' dat <- fixnames(dat, "query")
#'
#' # Define colors
#' map_gist(data=dat, color=c('#976AAE','#6B944D','#BD5945'))
#' map_gist(data=dat$gbif, color=c('#976AAE','#6B944D','#BD5945'))
#' map_gist(data=dat$ecoengine, color=c('#976AAE','#6B944D','#BD5945'))
#'
#' # Define colors and marker size
#' map_gist(data=dat, color=c('#976AAE','#6B944D','#BD5945'), size=c('small','medium','large'))
#'
#' # Define symbols
#' map_gist(data=dat, symbol=c('park','zoo','garden'))
#' }

map_gist <- function(data, description = "", public = TRUE, browse = TRUE, ...) {
  stopifnot(is(data, "occdatind") | is(data, "occdat"))
  data <- occ2df(data)
  datgeojson <- style_geojson(input = data, var = "name", ...)
  file <- tempfile(fileext = ".csv")
  write.csv(datgeojson, file)
  geofile <- togeojson2(file)
  gist_create(geofile, description = description, public = public, browse = browse)
}
