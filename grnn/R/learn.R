#' Learn
#' 
#' Create or update a General regression neural network.
#' 
#' @param set Data frame representing the training set. The first column is used to define the category of each observation (set \code{category.column} if it is not the case).
#' @param nn A General regression neural network with or without training.
#' @param variable.column The field number of the variable (1 by default).
#' @seealso \code{\link{grnn-package}}
#' @export
learn <- function(set, nn, variable.column=1) {
    if(missing(set)) stop("Set is missing!")
    if(missing(nn)) nn <- create.grnn()
    if(is.null(nn$set)) {
        nn$variable.column <- variable.column
        nn$set <- set
    } else {
        nn$set <- rbind(nn$set, set)
    }
    nn$Xa <- as.matrix(nn$set[,-nn$variable.column])
    nn$Ya <- as.matrix(nn$set[,nn$variable.column])
    nn$k <- length(nn$set[1,]) - 1
    nn$n <- length(nn$set[,1])
    return(nn)
}
