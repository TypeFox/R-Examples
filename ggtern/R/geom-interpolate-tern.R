#' Ternary Interpolation
#' 
#' This is the heavily requested geometry for interpolating between ternary values, results being
#' rendered using contours on a ternary mesh. 
#' 
#' @section Aesthetics: 
#' \Sexpr[results=rd,stage=build]{ggtern:::rd_aesthetics("geom", "InterpolateTern")}
#' @inheritParams geom_confidence_tern
#' @inheritParams ggplot2::geom_smooth
#' @inheritParams ggplot2::geom_point
#' @inheritParams ggplot2::geom_path
#' @examples 
#' data(Feldspar)
#' ggtern(Feldspar,aes(Ab,An,Or,value=T.C)) + 
#' stat_interpolate_tern(geom="polygon",
#'                      formula=value~x+y,
#'                      method=lm,n=100,
#'                      breaks=seq(0,1000,by=100),
#'                      aes(fill=..level..),expand=1) +
#'                      geom_point()
#' @name geom_interpolate_tern
#' @rdname geom_interpolate_tern
#' @export
geom_interpolate_tern <- function( mapping = NULL, data = NULL, stat = "InterpolateTern", position = "identity", 
                                   ...,
                                   lineend = "butt",linejoin = "round", linemitre = 1,
                                   na.rm = FALSE, show.legend = NA,inherit.aes = TRUE, 
                                   method='auto', formula=value~poly(x,y,degree=1)) {
  layer(
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomInterpolateTern,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(
      lineend     = lineend,
      linejoin    = linejoin,
      linemitre   = linemitre,
      na.rm       = na.rm,
      formula     = formula,
      method      = method,
      ...
    )
  )
}

#' @name geom_interpolate_tern
#' @rdname geom_interpolate_tern
#' @format NULL
#' @usage NULL
#' @export
GeomInterpolateTern <- ggproto("GeomInterpolateTern", 
                                GeomPath,
                                default_aes = aes(colour = "#3366FF", size = 0.5, linetype = 1, alpha = NA)
                                #draw_key = draw_key_path
)
