#' theme_white_f is a ggplot2 theme function that can be added to a ggplot2 object  
#' to eliminate axes, ticks and put white background
#' 
#' @param base_size base font size
#' @param base_family base font family
#' @export
theme_white_f <- function(base_size = 12, base_family = ""){
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      #  strip.background = element_blank(),
      plot.background = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}



#' Prepares the base of a map, starting from a data frame with longitude (long) and latitude (lat) coordinates, as a ggplot2 object
#' 
#' @param data data frame with longitude, latitude and group coordinates of a territory (polygons)
#' @param mapping Aesthetic mappings, as character strings (\code{link{ggplot2::aes_string}}). Defaults are \code{x = "long"}, \code{y = "lat"} and \code{group = "group"} (these work with every element of the list \code{shp_df} of ggswissmaps)
#' 
#'  
#' @examples 
#' data(shp_df)
#' maps2_(data = shp_df[[1]])
#' @export
maps2_ <- 
  function(data,
           mapping = 
             ggplot2::aes_string(x = "long",
                                 y = "lat",
                                 group = "group")){
    ggplot2::ggplot(data = data, mapping = mapping) +
    ggplot2::geom_path() +
    ggplot2::coord_equal() +
      ggplot2::theme(
        legend.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        #  strip.background = element_blank(),
        plot.background = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
}

# data(shp_df)
# maps2_(data = shp_df[[1]])
# maps2 <- lapply(shp_df, maps2_)