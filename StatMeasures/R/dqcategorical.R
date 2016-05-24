#' Data quality check of categorical variables
#' 
#'  Takes in a data, and returns summary of categorical variables
#'  @param data a data.frame or data.table
#'  @details
#'  While trying to understand a data, it is important to know the distribution of 
#'  categorical variables. \code{dqcategorical} produces an output which answers a 
#'  couple of questions regarding such variabes - how many distinct categories does 
#'  the variable have, what are the categories, what is the frequency
#'  of each of them and the percentage frequency.
#'  
#'  But first, it is critical to identify categorical variables in the data. They may be
#'  integer, numeric or character. All such variables should be converted to factor; one
#'  may use \code{factorise} function in this package to do this task easily.
#'  
#'  The function identifies all the factor variables and produces an output for each
#'  of them and returns a consolidated summary. It works for both 'data.frame' and 
#'  'data.table' but the output summary is a 'data.frame' only.
#'  @return a data.frame which contains the variable, category index, category, 
#'          category frequency and percentage frequency of all factor variables
#'  @author Akash Jain
#'  @seealso \code{\link{dqcontinuous}}, \code{\link{dqdate}}, \code{\link{contents}}
#'  @examples
#'  # A 'data.frame'
#' df <- data.frame(phone = c('IP', 'SN', 'HO', 'IP', 'SN', 'IP', 'HO', 'SN', 'IP', 'SN'),
#'                  colour = c('black', 'blue', 'green', 'blue', 'black', 'silver', 'black', 
#'                  'white', 'black', 'green'))
#'                  
#' # Factorise categorical variables
#' df <- factorise(data = df, colNames = c('phone', 'colour'))
#' 
#' # Generate a data quality report of continuous variables
#' summaryCategorical <- dqcategorical(data = df)
#'  @export
dqcategorical <- function(data) {
  if(class(data)[1] != 'data.frame' && class(data)[1] != 'data.table') {
    stop('Invalid input: data should be either data.frame or data.table')
  } else {
    varNames <- names(data)
    classVar <- sapply(data, class)
    catVars <- varNames[classVar == 'factor']
    len <- length(catVars)
    if(len > 0) {
      if(class(data)[1] == 'data.frame') {
        dataCatVars <- data[catVars]
      } else if (class(data)[1] == 'data.table'){
        dataCatVars <- data[, catVars, with = FALSE]
      }
      catVarsFreqList <- lapply(dataCatVars, function(var) data.frame(table(var, useNA = 'ifany')))
      catVarsFreqSummary <- do.call(rbind, catVarsFreqList)
      variable <- sapply(strsplit(row.names(catVarsFreqSummary), '.', fixed = TRUE), function(vec) vec[1])
      categorynum <- sapply(strsplit(row.names(catVarsFreqSummary), '.', fixed = TRUE), function(vec) vec[2])    
      catVarsFreqSummary <- data.frame(variable,
                                       categorynum,
                                       catVarsFreqSummary, 
                                       row.names = NULL)
      catVarsFreqSummary[, 'percent'] <- round((catVarsFreqSummary[, 'Freq']/nrow(data))*100, digits = 2)
      names(catVarsFreqSummary) <- c('variable', 'categorynum', 'category', 'frequency', 'percentage')
      return(catVarsFreqSummary)
    } else {
      stop('There are no variables of class factor in the data')
    }
  }
}