#' interval.check
#' 
#' Check which interval a number belongs to
#' 
#' This function takes in a data.frame with a specified column and compares that to a vector of times
#' 
#' @author Jared P. Lander
#' @aliases interval.check
#' @export interval.check
#' @param data data.frame
#' @param input character name of column we wish to compare
#' @param times vector in ascending order where the differences between sequential elements are the intervals
#' @param fun character containing compator
#' @return Vector indicating which element of \code{times} that row belongs to.  If the row is beyond any element NA is in it's spot.
#' @examples
#' 
#' head(cars)
#' interval.check(cars, input="speed", times=seq(min(cars$speed), max(cars$speed), length=10))
interval.check <- function(data, input="Stop", times, fun="<=")
{
    # do an outer product seeing which of input is meets the requirement in realtion to each of the times
    equalityMat <- outer(data[, input], times, FUN=fun)
    
    ## now we are going to see the first of the columns in each row to hold TRUE
    ## the best way is to do a cumsum, by row, of the negation of each cell
    equalityMat <- t(apply(!equalityMat, 1, cumsum))    # transpose it get it back in the right shape
    
    # now take the max of each row, this will tell you the first column that it worked for
    # add one to account that the first column is 0
    indices <- apply(equalityMat, 1, max) + 1
    
    # now determine which upper bound it belongs to, NA means it didn't happen within alloted time
    belongsTo <- times[indices]
    
    return(belongsTo)
}
