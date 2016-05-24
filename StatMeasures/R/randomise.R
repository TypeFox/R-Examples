#' Order the rows of a data randomly
#' 
#'  Takes in data and seed, and returns the data with randomly ordered observations
#'  @param data a matrix, data.frame or data.table
#'  @param seed an integer value
#'  @details
#'  Some of the modeling algorithms pick top p percent of the observations for
#'  training the model, which could lead to skewed predictions. This function
#'  solves that problem by randomly ordering the observations so that the response
#'  variable has more or less the same distribution even if the algorithms don't pick
#'  training observations randomly.
#'  @return \code{data} of same class as input with randomly ordered observations
#'  @author Akash Jain
#'  @seealso \code{\link{factorise}}, \code{\link{rmdupkey}}, \code{\link{rmdupobs}}
#'  @examples
#'  # A 'data.frame'
#' df <- data.frame(x = c(1, 2, 3, 4, 5), y = c('a', 'b', 'c', 'd', 'e'))
#' 
#' # Change the order of the observations randomly
#' dfRan <- randomise(data = df)
#' dfRan <- randomise(data = df, seed = 150)
#'  @export
randomise <- function(data, seed = NULL) {
  if(class(data)[1] != 'matrix' && class(data)[1] != 'data.frame' && class(data)[1] != 'data.table') {
    stop('Invalid input: data should be either matrix or data.frame or data.table')
  } else {
    if(!is.null(seed)) set.seed(seed)
    n <- nrow(data)
    data <- data[sample(1:n, size = n), ]
    rownames(data) <- NULL
    return(data)    
  }
}