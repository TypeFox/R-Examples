#' Basic summary of the data
#' 
#'  Takes in a data and returns summary of the data
#'  @param data a data.frame or data.table
#'  @details
#'  This function helps when one wants to get a quick snapshot of the data such as class,
#'  distinct values, missing values and sample value of the variables.
#'  
#'  It works for both 'data.frame' and 'data.table' but the output will be a 'data.frame' only.
#'  @return a data.frame that contains variable, class, distinct values,
#'          missing values, percentage of missing value and sample value
#'  @author Akash Jain
#'  @seealso \code{\link{dqcontinuous}}, \code{\link{dqcategorical}}, \code{\link{dqdate}}
#'  @examples
#'  # A data frame
#' df <- data.frame(x = c(1, 2, 3, 4, NA),
#'                  y = c('a', 'b', 'c', 'd', 'e'),
#'                  z = c(1, 1, 0, 0, 1))
#'
#' # Summary of the data
#' dfContents <- contents(data = df)
#'  @export
contents <- function(data) {
  if(class(data)[1] != 'data.frame' && class(data)[1] != 'data.table') {
    stop('Invalid input: data should be either data.frame or data.table')
  } else {
    variable <- names(data)
    class <- sapply(data, class)
    distinctValues <- sapply(data, function(var) length(unique(var)))
    missingValues <- sapply(data, function(var) sum(is.na(var)))
    perMissingValues <- sapply(missingValues, function(value) round(value/nrow(data), digits = 4)*100)
    sampleValue <- sapply(data, function(var) sample(var, size = 1))
    contents <- data.frame(variable, 
                           class, 
                           distinctValues, 
                           missingValues, 
                           perMissingValues, 
                           sampleValue,
                           row.names = NULL)
    return(contents)    
  }
}