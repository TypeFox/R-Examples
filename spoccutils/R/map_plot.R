#' Use base R plotting method to visualize spocc data.
#'
#' @importFrom rworldmap getMap
#' @export
#' @param x Input, an object of class \code{occdat}
#' @param ... Further args to \code{\link{points}}
#' @examples \dontrun{
#' library("spocc")
#' (out <- occ(query='Accipiter striatus', from='gbif', limit=100, has_coords=TRUE))
#' map_plot(out)
#' }
map_plot <- function(x, ...) {
  df <- occ2df(x)
  df <- df[complete.cases(df),]
  df <- df[df$longitude != 0, ]
  sp::coordinates(df) <- ~longitude + latitude
  sp::proj4string(df) <- sp::CRS("+init=epsg:4326")
  sp::plot(rworldmap::getMap())
  points(df, col = "red", ...)
}
