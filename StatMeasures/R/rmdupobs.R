#' Remove duplicate observations from data
#' 
#'  Takes in a data, and returns it with duplicate observations removed
#'  @param data a data.frame or data.table
#'  @details
#'  Duplicate observations are redundant and they need to be removed from the data.
#'  \code{rmdupobs} does just that; it removes the duplicated observations (the ones
#'  in which value of every variable is duplicated) and returns the data with only 
#'  unique observations.
#'  
#'  It works for both 'data.frame' and 'data.table' and returns the \code{data} with 
#'  same class as that of input.
#'  @return a \code{data} of same class as input with only unique observations
#'  @author Akash Jain
#'  @seealso \code{\link{randomise}}, \code{\link{rmdupkey}}, \code{\link{factorise}}
#'  @examples
#'  # A 'data.frame'
#' df <- data.frame(x = c(1, 2, 5, 1), y = c(3, 3, 1, 3))
#'
#' # Remove duplicate observations from data
#' dfUnq <- rmdupobs(data = df)
#'  @export
rmdupobs <- function(data) {
  if(class(data)[1] != 'data.frame' && class(data)[1] != 'data.table') {
    stop('Invalid input: data should be either data.frame or data.table')
  } else if(class(data)[1] == 'data.frame') {
    unqData <- data[!duplicated(data), ]
  } else if(class(data)[1] == 'data.table') {
    unqData <- unique(data)
  }
  diff <- nrow(data) - nrow(unqData)
  return(unqData)
}