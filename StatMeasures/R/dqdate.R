#' Data quality check of date variables
#' 
#'  Takes in a data, and returns summary of date variables
#'  @param data a data.frame or data.table
#'  @details
#'  \code{dqdate} produces summary of all date variables in the data. The function
#'  identifies all variables as date if they are of class 'Date' or 'IDate'.
#'  
#'  Generally the dates are imported in R as character. They must be converted to
#'  an appropriate date format and then the function should be used.
#'  
#'  The summary includes variable, non-missing values, missing values, minimum and
#'  maximum of the date variabes. Input data can be a 'data.frame' or 'data.table' but the
#'  output summary will be a 'data.frame' only.
#'  @return a data.frame which contains the variable, non-missing values, missing values, 
#'          minimum and maximum of all date variables
#'  @author Akash Jain
#'  @seealso \code{\link{dqcontinuous}}, \code{\link{dqcategorical}}, \code{\link{contents}}
#'  @examples
#'  # A 'data.frame'
#' df <- data.frame(date = c('2012-11-21', '2015-1-4', '1996-4-30', '1979-9-23', '2005-5-13'),
#'                  temperature = c(26, 32, 35, 7, 14))
#'                  
#' # Convert character date to date format
#' df[, 'date'] <- as.Date(df[, 'date'])
#' 
#' # Generate a data quality report of date variables
#' summaryDate <- dqdate(data = df)
#'  @export
dqdate <- function(data) {
  if(class(data)[1] != 'data.frame' && class(data)[1] != 'data.table') {
    stop('Invalid input: data should be either data.frame or data.table')
  } else {
    varNames <- names(data)
    classVar <- sapply(data, function(var) class(var)[1])
    dateVars <- varNames[classVar == 'Date' | classVar == c('IDate')]
    len <- length(dateVars)
    if(len > 0) {
      dataDateVar <- data[dateVars]
      variable <- names(dataDateVar)
      nonMissingValues <- sapply(dataDateVar, function(var) sum(!is.na(var)))
      missingValues <- sapply(dataDateVar, function(var) sum(is.na(var)))
      maximum <- lapply(dataDateVar, function(var) max(var, na.rm=T))
      minimum <- lapply(dataDateVar, function(var) min(var, na.rm=T))
      dateVarSummary <- data.frame(variable, 
                                   nonMissingValues, 
                                   missingValues, 
                                   minimum, 
                                   maximum, 
                                   row.names=NULL)
      names(dateVarSummary)[c(4,5)] <- c('minimum', 'maximum')
      return(dateVarSummary)
    } else {
      stop('There are no variables of class Date or IDate in the data')
    }
  }
}