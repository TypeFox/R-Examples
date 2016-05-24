#' Convert adpcr to ppp
#' 
#' Converts \code{\linkS4class{adpcr}} object to the list of
#' \code{\link[spatstat]{ppp.object}}s allowing spatial analysis.
#' 
#' @details 
#' Each array is independently converted by \code{\link[spatstat]{ppp}}
#' function. \code{marks} attached to each point represent values contained by
#' the \code{\linkS4class{adpcr}} object.
#' 
#' @param input Object of the \code{\linkS4class{adpcr}} class containing data
#' from one or more panels.
#' @param marks If \code{TRUE}, marks values for non-empty partitions.
#' @param plot If \code{TRUE}, array is plotted.
#' @return A list containing objects with class
#' \code{\link[spatstat]{ppp.object}} with the length equal to the number of
#' arrays (minimum 1).
#' @author Michal Burdukiewcz, Stefan Roediger.
#' @seealso \code{\link[spatstat]{ppp.object}}, \code{\link[spatstat]{ppp}}.
#' @keywords manip panel
#' @export adpcr2ppp
#' @examples
#' 
#' many_panels <- sim_adpcr(m = 400, n = 765, times = 1000, pos_sums = FALSE, 
#'                    n_panels = 5)
#' 
#' # Convert all arrays to ppp objects
#' adpcr2ppp(many_panels)
#' 
#' # Convert all arrays to ppp objects and get third plate
#' third_plate <- adpcr2ppp(many_panels)[[3]]
#' 
#' # Convert only third plate to ppp object
#' third_plate2 <- adpcr2ppp(extract_dpcr(many_panels, 3))
#' 
#' # Check the class of a new object
#' class(third_plate2)
#' 
#' # It's a list with the length 1. The third plate is a first element on this 
#' #list
#' class(third_plate2[[1]])
#' 
#' 
adpcr2ppp <- function(input, marks = TRUE, plot = FALSE) {
  arrays <- adpcr2panel(input)
  
  lapply(arrays, function(single_array)
    create_ppp(data_vector = single_array, nx_a = ncol(single_array), ny_a = nrow(single_array),
               marks = marks, plot = plot))
}


create_ppp <- function(data_vector, nx_a, ny_a, plot, marks) {
  #strange syntax, because spatstat use different localizations
  #than dpcR.
  data_points <- which(matrix(data_vector, ncol = nx_a, nrow = ny_a) > 0,
                       arr.ind = TRUE)
  data_points[, "row"] <- ny_a - data_points[, "row"] + 1
  
  if (plot)
    plot(ppp(data_points[, 2], data_points[, 1], c(1, nx_a), c(1, ny_a)))
  #check if marks are properly assigned
  if (marks) {
    ppp(data_points[, 2], data_points[, 1], c(1, nx_a), c(1, ny_a), 
        marks = data_vector[data_vector != 0])
  } else {
    ppp(data_points[, 2], data_points[, 1], c(1, nx_a), c(1, ny_a))
  }
}