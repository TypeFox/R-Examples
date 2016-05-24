#' shift.column
#' 
#' Shift a column of data
#' 
#' Shifts a column of data up or down a certain number of rows
#' 
#' @author Jared P. Lander
#' @aliases shift.column
#' @export shift.column
#' @return \code{\link{data.frame}} with the specified columns shifted.
#' @param data \code{\link{data.frame}}
#' @param columns Character vector specifying which columns to shift.
#' @param newNames Character vector specifying new names for the columns that will be created by the shift.  Must be same length as \code{columns}.
#' @param len Integer specifying how many rows to shift the data.
#' @param up logical indicating if rows should be shifted up or down.
#' @examples
#' 
#' myData <- data.frame(Upper=LETTERS, lower=letters)
#' shift.column(data=myData, columns="lower")
#' shift.column(data=myData, columns="lower", len=2)
#' 
shift.column <- function(data, columns, newNames=sprintf("%s.Shifted", columns), len=1L, up=TRUE)
{
    if(length(columns) != length(newNames))
    {
        stop("columns and newNames must be the same length")
    }
    
    # get the rows to keep based on how much to shift it by and weather to shift up or down
    rowsToKeep <- seq(from=1 + len*up, length.out=NROW(data) - len)
    
    # for the original dat ait needs to be shifted the other way
    dataRowsToKeep <- seq(from=1 + len*!up, length.out=NROW(data) - len)
    
    #create a df of the shifted rows
    shiftedDF <- data[rowsToKeep, columns]
    
    # give the right names to these new columns
    names(shiftedDF) <- newNames
    
    # data names
    dataNames <- names(data)
    
    # get rid of excess rows in data
    data <- data[dataRowsToKeep, ]
    
    # tack shifted data onto the end of the original (and cutoff) data
    data <- cbind(data, shiftedDF)
    names(data) <- c(dataNames, newNames)
    
    # return the data
    return(data)
}
