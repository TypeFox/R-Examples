#' Render Borders on Top
#' 
#' Convenience functions to render the axis border lines on top (or bottom) of the other layers. 
#' By default the borders are rendered in the background (bottom)
#' @author Nicholas Hamilton
#' @rdname theme_bordersontop
#' @export
theme_bordersontop = function(){
  theme(tern.axis.line.ontop=TRUE)
}

#' @rdname theme_bordersontop
#' @export
theme_bordersonbottom = function(){
  theme(tern.axis.line.ontop=FALSE)
}