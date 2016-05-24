#' Learn
#' 
#' Create or update a Probabilist neural network.
#' 
#' The function \code{learn} aims to create a new Probabilist neural network with a training set, or update the training set of an already trained Probabilist neural network. It sets the parameters \code{model}, \code{set}, \code{category.column}, \code{categories}, \code{k} and \code{n} of the neural network.
#' 
#' @param set Data frame representing the training set. The first column is used to define the category of each observation (set \code{category.column} if it is not the case).
#' @param nn A Probabilistic neural network with or without training.
#' @param category.column The field number of the category (1 by default).
#' 
#' @seealso \code{\link{pnn-package}}, \code{\link{smooth}}, \code{\link{perf}}, \code{\link{guess}}, \code{\link{norms}}
#' 
#' @export
#' 
#' @return A trained Probabilist neural network.
#' 
#' @examples
#' library(pnn)
#' data(norms)
#' pnn <- learn(norms)
#' pnn$model
#' pnn$set[1:10,]
#' pnn$category.column
#' pnn$categories
#' pnn$k
#' pnn$n
learn <- function(set, nn, category.column=1) {
    if(missing(set)) stop("Set is missing!")
    if(missing(nn)) nn <- create.pnn()
    if(is.null(nn$set)) {
        nn$category.column <- category.column
        nn$set <- set
    } else {
        nn$set <- rbind(nn$set, set)
    }
    nn$set[,nn$category.column] <- factor(nn$set[,nn$category.column])
    nn$categories <- levels(nn$set[,nn$category.column])
    nn$k <- length(nn$set[1,]) - 1
    nn$n <- length(nn$set[,1])
    # Scale
    return(nn)
}
