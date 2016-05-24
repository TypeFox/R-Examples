




#' Offers various swiss maps as ggplot2 objects.
#'
#' Offers various swiss maps as ggplot2 objects and gives the possibility
#' to add layers of data on the maps. Data are publicly available from the
#' swiss federal statistical office.
#'
#' There are only two objects in package \pkg{ggswissmaps}:
#' \code{theme_white} (for white background of ggplot2 objects) and \code{maps} (a list of 8 swiss maps).
#' 
#' 
#' @docType package
#' @name ggswissmaps
#' 
#' @import ggplot2
NULL


#' theme_white is a ggplot2 theme object and can be added to a ggplot2 object  
#' to eliminate axes, ticks and put white background
theme_white <- ggplot2::theme(
  axis.title=element_blank(),
  axis.text=element_blank(),
  axis.ticks=element_blank(),
  axis.line=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.margin=element_blank(),
  panel.grid=element_blank()
)
#' @name theme_white
NULL


#' A list with 16 data frames of swiss territory boundaries, at various levels.
#'
#' Every element of the list is a data frame, which can be used with ggplot2.
#'
#' Columns are not all the same across data frames, but usually they have the following in common:
#' \itemize{
#'   \item long. Longitude coordinate (x)
#'   \item lat. Latitude coordinate (y)
#'   \item group. A factor to be used to plot the polygons correctly (with ggplot2)
#' }
#'
#' @format A list with 16 data frames with swiss territory boundaries (at various levels).
#' @source \url{http://www.bfs.admin.ch/bfs/portal/fr/index/dienstleistungen/geostat/datenbeschreibung.html}
#' @examples
#' data(shp_df)
#' class(shp_df)
#' length(shp_df)
#' names(shp_df)
#' str(shp_df[["g1k15"]])
#' @name shp_df
NULL

#' A list with 16 ggplot2 objects of swiss territory boundaries, at various levels.
#'
#' Every element of the list is a ggplot2 graphic, corresponding to an element of \code{shp_df}.
#'
#' Columns are not all the same across data frames, but usually they have the following in common:
#' \itemize{
#'   \item long. Longitude coordinate (x)
#'   \item lat. Latitude coordinate (y)
#'   \item group. A factor to be used to plot the polygons correctly (with ggplot2)
#' }
#'
#' @format A list with 16 data frames with swiss territory boundaries (at various levels).
#' @source \url{http://www.bfs.admin.ch/bfs/portal/fr/index/dienstleistungen/geostat/datenbeschreibung.html}
#' @examples
#' class(maps2)
#' length(maps2)
#' names(maps2)
#' # str(maps2[["g1k15"]])
#' 
#' # By name
#' maps2[["g1k15"]]
#' 
#' # By index
#' maps2[[5]]
#' @export
#' @name maps2
data("shp_df", envir = parent.env(environment()))
maps2 <- lapply(shp_df, maps2_)

# Sostituito le due righe qui sotto (0.0.7) (e tolto maps2.rda dalla cartella data/)
# #' @name maps2
# NULL

# Con queste 4 (0.0.8)
# #' @name maps2
# #' @export
# data("shp_df")
# maps2 <- lapply(shp_df, maps2_)


# extract_data <- function(gg){
#   l <- ggplot2::ggplot_build(gg)
#   l$data[[1]]
# }
# 
# x <- extract_data(maps[[1]])
# str(x)
# 
# x <- extract_data(maps[["cantoni"]])
# str(x)
# 
# pcant <- ggplot(data = x, mapping = aes(x=x, y=y, group=group, fill=factor(group))) + geom_polygon() + coord_equal() + theme_white
# pcant
# pcant + theme(legend.position="none")