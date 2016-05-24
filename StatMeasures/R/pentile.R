#' Create pentiles of a variable
#' 
#'  Takes in a vector, and returns a vector of pentiles
#'  @param vector an integer or numeric vector
#'  @param decreasing a logical input, which if set to \code{FALSE} puts smallest values in
#'         pentile 1 and if set to \code{TRUE} puts smallest values in pentile 5; \code{FALSE}
#'         is default
#'  @details
#'  \code{pentile} is a convinient function to get integer pentiles of an integer or 
#'  numeric vector. By default, the smallest values are placed in the smallest pentile.
#'  
#'  Sometimes one may want to put smallest values in the biggest pentile, and for that 
#'  the user can set the \code{decreasing} argument to \code{TRUE}; by default it is
#'  \code{FALSE}.
#'  @return an integer vector of pentile values
#'  @author Akash Jain
#'  @seealso \code{\link{decile}}, \code{\link{outliers}}, \code{\link{imputemiss}}
#'  @examples
#'  # Scores vector
#' scores <- c(1, 4, 7, 10, 15, 21, 25, 27, 32, 35, 
#'             49, 60, 75, 23, 45, 86, 26, 38, 34, 67)
#' 
#' # Create pentiles based on the values of the vector
#' pentileScores <- pentile(vector = scores)
#' pentileScores <- pentile(vector = scores, decreasing = TRUE)
#'  @export
pentile <- function(vector, decreasing = FALSE) {
  if(class(vector) != 'integer' && class(vector) != 'numeric') {
    stop('Invalid input: vector should be either integer or numeric')    
  } else if (class(decreasing) != 'logical' | length(decreasing) != 1) {
    stop('Invalid input: decreasing should be a logical vector of length 1')
  } else {
    pentile <- as.integer(as.character(cut(vector, 
                                           breaks = quantile(vector, probs = seq(0, 1, by = 0.2), na.rm = TRUE), 
                                           labels = if(decreasing == FALSE) {1:5} else {5:1}, 
                                           include.lowest = TRUE)))
    return(pentile)
  }
}