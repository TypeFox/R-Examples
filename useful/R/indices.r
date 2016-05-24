#' @title indexToPosition
#' @description Given a long matrix index convert to row and column positions
#' @details Using \code{\link{which}} on a matrix returns a number that iterates down rows then across columns.  This function returns the (row, column) position of that index.
#' @author Jared P. Lander
#' @aliases indexToPosition
#' @export indexToPosition
#' @param x Position of indices
#' @param nrow The number of rows in the matrix
#' @return A Matrix with row and column columns and a row for each value of \code{x}
#' @examples 
#' indexToPosition(3, 2)
#' indexToPosition(c(1, 4, 5, 7, 9), 3)
#' indexToPosition(1:16, 4)
#' indexToPosition(c(1, 3, 5, 6, 8, 10, 11, 13, 15), 5)
#' 
indexToPosition <- function(x, nrow=1)
{
    # first find the column
    theCol <- ceiling(x / nrow)
    # get the row
    theRow <- x - (theCol - 1)*nrow
    
    return(cbind(row=theRow, col=theCol))
}

#' @title positionToIndex
#' @description Given row and column positions calculate the index.
#' @details With row and column positions this computes the index, starting at (1,1) working down rows then across columns.
#' @author Jared P. Lander
#' @aliases positionToIndex
#' @export positionToIndex
#' @param row Vector specifying row positions
#' @param col Vector specifying column positions
#' @param nrow The number of rows in the matrix
#' @return A vector of indices
#' @examples 
#' positionToIndex(1, 2, 2)
#' positionToIndex(row=c(1, 1, 2, 1, 3), col=c(1, 2, 2, 3, 3), nrow=3)
#' positionToIndex(rep(1:4, 4), rep(1:4, each=4), nrow=4)
#' positionToIndex(rep(c(1, 3, 5), 3), rep(1:3, each=3), nrow=5)
#' 
positionToIndex <- function(row, col, nrow=max(row))
{
    (col - 1)*nrow + row
}