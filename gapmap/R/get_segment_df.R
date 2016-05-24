#'Make a data.frame object
#'
#'This function just make a data.frame based on 4 input parameters
#'
#' @param x0 x coordinate of point 1
#' @param y0 y coordinate of point 1
#' @param x1 x coordinate of point 2
#' @param y1 y coordinate of point 2
#' @export get_segment_df
#' @aliases get_segment_df
#' @return A data.frame
#' @keywords internal
#' 



#return a data frame object
get_segment_df <- function(x0, y0, x1, y1) {
  data.frame(x0, y0, x1, y1)
}