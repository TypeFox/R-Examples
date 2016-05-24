#'@title Create an Index of lowest values of each column for a matrix
#'
#'@description This function creates a vector filled with the row position
#'of the lowest value within each column of a matrix.
#'
#'@usage minindc(x)
#'
#'@param x A matrix that will be be used to create the index of lowest values.
#'
#'@return pos A vector containig an index of the row position of the lowest
#' value within each column of a matrix.
#'
#'

minindc <- function(x) {
    
    ncols <- dim(x)[2]
    
    nrows <- dim(x)[1]
    
    pos <- matrix(0, ncols, 1)
    
    for (i in 1:ncols) {
        
        pos[i, ] <- which.min(x)
    }
    return(pos)
} 
