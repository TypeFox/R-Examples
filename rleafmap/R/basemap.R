#'Define a Tile Basemap Layer
#'
#'Define a new basemap layer from a tile server.
#'
#'@param URL a character string giving the tile server url or a the name of a pre-configured server.
#'@param name a character string to name the layer.
#'@param alpha a numeric value in \eqn{[0, 1]} setting the layer opacity.
#'@param minZoom,maxZoom numeric values setting the minimum and maximum zoom level.
#'@param tileSize a numeric value setting tile size (width and height in pixels, assuming tiles are square).
#'@param tms logical. If \code{TRUE}, inverses Y axis numbering for tiles (for TMS services)
#'
#'@details \code{URL} should have the form \code{'http://{s}.somedomain.com/somepath/{z}/{x}/{y}.png'}
#'with \code{{s}} a facultative subdomain, \code{{z}} the zoom level and \code{{x}}, \code{{y}} the coordinates.
#'\pkg{rleafmap} comes with a list of pre-configured servers. Names of these servers are returned
#'by the function \code{bmSource}.
#'
#'@return An object of class \code{basemap} which can be directly used in \code{\link{writeMap}}.
#'@seealso \code{\link{spLayer}} to define data layers.
#'@examples \dontrun{
#'  #A simple map with two nice basemaps.
#'  bm1 <- basemap("mapquest.map")
#'  bm2 <- basemap("stamen.watercolor")
#'  writeMap(bm1, bm2)
#'}
#'@export
basemap <- function(URL, name=NULL, alpha=1, minZoom=0, maxZoom=18, tileSize=256, tms=FALSE){
  if(is.vector(URL) && length(URL)==1){
    URL <- as.character(URL)
    URL <- bmServer(URL)
    URL.cr <- bmCredit(URL)
  }else{
    stop("URL must be a single value character vector")
  }
  if(!is.numeric(alpha)){
    stop("alpha must be numeric")
  }else{
    if (alpha<0 || alpha>1){
      stop("alpha must be comprise between 0 and 1")
    }
  }
  if(!is.numeric(minZoom))
    stop("minZoom must be numeric")
  if(!is.numeric(maxZoom))
    stop("maxZoom must be numeric")
  if(!is.logical(tms))
    stop("tms must be set on TRUE or FALSE")  
  
  res <- list(name, URL, URL.cr, alpha, minZoom, maxZoom, tileSize, tolower(as.character(tms)))
  names(res) <- c("name", "URL", "URL.cr", "alpha", "minZoom", "maxZoom", "tileSize", "tms")
  class(res) <- "basemap"
  return(res) 
}
