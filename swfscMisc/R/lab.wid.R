#' @title Label Width
#' @description Calculate width of labels for plots.
#'
#' @param labels vector of labels to be used on plots
#' 
#' @importFrom graphics par strwidth
#' @export
#' 
lab.wid <- function(labels) {
  inches_to_lines <- (par("mar") / par("mai"))[1]
  max(strwidth(labels, units = "inches")) * inches_to_lines
}
