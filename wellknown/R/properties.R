#' Add properties to a geojson object
#'
#' @export
#'
#' @param x Input
#' @param style List of color, fillColor, etc., or NULL
#' @param popup Popup text, or NULL
#' @examples
#' str <- "POINT (-116.4000000000000057 45.2000000000000028)"
#' x <- wkt2geojson(str)
#' properties(x, style=list(color = "red"))
properties <- function(x, style = NULL, popup = NULL){
  if(is.null(style) && is.null(popup)) {
    stop("You must supply a list of named options to either style, popup, or both", call. = FALSE)
  } else {
    if(is(style, "list") || is(popup, "list")){
      if(length(style) == 0 && length(popup) == 0)
        stop("At least one of style or popup needs a non-empty list", call. = FALSE)
    }
    modifyList(x, list(properties = list(style = style, popup = popup)))
  }
}
