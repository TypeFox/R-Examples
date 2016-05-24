#' ggplot2 visualization of species occurences
#'
#' @importFrom grid unit
#' @export
#' @param x Input, object of class \code{occdat}
#' @param map (character) One of world, world2, state, usa, county, france, italy, or nz
#' @param point_color Default color of your points
#' @examples \dontrun{
#' library("spocc")
#' dat <- occ(query = 'Lynx rufus californicus', from = 'ecoengine', limit=100)
#' map_ggplot(dat)
#' map_ggplot(dat, "usa")
#' map_ggplot(dat, "county")
#'}
map_ggplot <- function(x, map = "world", point_color = "#86161f") {
  latitude <- longitude <- lat <- long <- group <- NA
  dt <- occ2df(x)
  dt <- dt[complete.cases(dt), ]
  wmap <- map_data(map)
  ggplot(dt, aes(longitude, latitude)) +
    geom_point(color = point_color, size = 3) +
    geom_polygon(aes(long, lat, group = group), fill = NA, colour = "black", data = wmap) +
    sutils_blank_theme()
}

sutils_blank_theme <- function(){
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        plot.margin = rep(unit(0, "null"), 4))
}
