#' Split modeling data into test and train set
#' 
#'  Takes in data, fraction (for train set) and seed, and returns train and test set
#'  @param data a matrix, data.frame or data.table
#'  @param fraction proportion of observations that should go in the train set
#'  @param seed an integer value
#'  @details
#'  An essential task before doing modeling is to split the modeling data into
#'  train and test sets. \code{splitdata} is built for this task and returns a list
#'  with train and test sets, which can be picked using the code given in example.
#'  
#'  \code{fraction} corresponds to the train dataset, while the rest of the
#'  observations go to the test dataset. If the user wants to generate the same
#'  test and train dataset everytime, he should specify a \code{seed} value.
#'  @return a list with two elements: train and test set
#'  @author Akash Jain
#'  @seealso \code{\link{actvspred}}, \code{\link{mape}}, \code{\link{accuracy}}, 
#'           \code{\link{auc}}, \code{\link{iv}}, \code{\link{ks}}
#'  @examples
#'  # A 'data.frame'
#' df <- data.frame(x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), 
#'                  y = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'),
#'                  z = c(1, 1, 0, 0, 1, 0, 0, 1, 1, 0))
#' 
#' # Split data into train (70%) and test (30%)
#' ltData <- splitdata(data = df, fraction = 0.7, seed = 123)
#' trainData <- ltData$train
#' testData <- ltData$test
#'  @export
splitdata <- function(data, fraction, seed = NULL) {
  if(class(data)[1] != 'matrix' && class(data)[1] != 'data.frame' && class(data)[1] != 'data.table') {
    stop('Invalid input: data should be either matrix or data.frame or data.table')
  } else if(fraction < 0 || fraction > 1) {
    stop('Invalid input: fraction should be in the range 0 and 1')
  } else {
    if(!is.null(seed)) set.seed(seed)
    index <- 1:nrow(data)
    trainIndex <- sample(index, size = floor(fraction * nrow(data)))
    train <- data[trainIndex, ]
    test <- data[-trainIndex, ]
    rownames(train) <- NULL
    rownames(test) <- NULL
    return(list(train = train, test = test))    
  }
}