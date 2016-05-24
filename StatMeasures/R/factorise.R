#' Change the class of variables to factor
#' 
#'  Takes in data and colNames, and returns the data with all variables mentioned in
#'  colNames converted to factor
#'  @param data a data.frame or data.table
#'  @param colNames a character vector of variable names to be converted to factor
#'  @details
#'  We often face the task of converting a bunch of variables to factor. This
#'  function is particularly useful in such a situation. Just specify the \code{data}
#'  and variable names and all of them will be converted to factor.
#'  
#'  It works for both 'data.frame' and 'data.table', and the output is data of the same
#'  class as that of input.
#'  @return \code{data} of same class as input with specified variables converted to factor
#'  @author Akash Jain
#'  @seealso \code{\link{randomise}}, \code{\link{rmdupkey}}, \code{\link{rmdupobs}}
#'  @examples
#'  # A 'data.frame'
#' df <- data.frame(x = c(1, 2, 3, 4, 5), 
#'                  y = c('a', 'b', 'c', 'd', 'e'),
#'                  z = c(1, 1, 0, 0, 1))
#' 
#' # Change the class of variables y and z to factors
#' dfFac <- factorise(data = df, colNames = c('y', 'z'))
#'  @export
#'  @importFrom data.table .SD
factorise <- function(data, colNames) {
  if(class(colNames) != 'character') {
    stop('Invalid input: colNames should be a character vector')    
  } else if(class(data)[1] != 'data.frame' && class(data)[1] != 'data.table') {
    stop('Invalid input: data should be either data.frame or data.table')
  } else if(sum(colNames %in% names(data)) != length(colNames)) {
    stop(paste('Column', paste(colNames[!colNames %in% names(data)], collapse = ', '), 'not present in the data'))
  } else {
    id <- which(names(data) %in% colNames)
    if(class(data)[1] == 'data.frame'){
      data[, id] <- lapply(data.frame(data[, id]), as.factor)
    } else if(class(data)[1] == 'data.table'){
      data[, id] <- data[, lapply(.SD, as.factor), .SDcols = id]
    }
    return(data)
  }
}