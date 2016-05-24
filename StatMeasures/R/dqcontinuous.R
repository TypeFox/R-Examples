#' Data quality check of continuous variables
#' 
#'  Takes in a data, and returns summary of continuous variables
#'  @param data a data.frame or data.table
#'  @details
#'  It is of utmost importance to know the distribution of continuous variables in the
#'  data. \code{dqcontinuous} produces an output which tells - continuous variable,
#'  non-missing values, missing values, percentage missing, minumum, average, maximum,
#'  standard deviation, variance, common percentiles from 1 to 99, and number of outliers
#'  for each continuous variable.
#'  
#'  The function tags all integer and numeric variables as continuous, and produces output
#'  for them; if you think there are some variables which are integer or numeric in the
#'  data but they don't represent a continuous variable, change their type to an 
#'  appropriate class.
#'  
#'  \code{dqcontinuous} uses the same criteria to identify outliers as the one used for
#'  box plots. All values that are greater than 75th percentile value + 1.5 times the 
#'  inter quartile range or lesser than 25th percentile value - 1.5 times the inter
#'  quartile range, are tagged as outliers.
#'  
#'  This function works for both 'data.frame and 'data.table' but returns a 'data.frame' only.
#'  @return a data.frame which contains the non-missing values, missing values,
#'          percentage of missing values, mimimum, mean, maximum, standard deviation,
#'          variance, percentiles and count of outliers of all integer and 
#'          numeric variables
#'  @author Akash Jain
#'  @seealso \code{\link{dqcategorical}}, \code{\link{dqdate}}, \code{\link{contents}}
#'  @examples
#'  # A 'data.frame'
#' df <- data.frame(x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#'                  y = c(22, NA, 66, 12, 78, 34, 590, 97, 56, 37))
#' 
#' # Generate a data quality report of continuous variables
#' summaryContinuous <- dqcontinuous(data = df)
#'  @export
dqcontinuous <- function(data) {
  if(class(data)[1] != 'data.frame' && class(data)[1] != 'data.table') {
    stop('Invalid input: data should be either data.frame or data.table')
  } else {
    numoutliers <- function(vector) {
      if(class(vector) != 'integer' && class(vector) != 'numeric') {
        stop('Invalid input: vector should be either integer or numeric')
      } else {
        p25 <- quantile(vector, c(0.25), na.rm = TRUE)
        p75 <- quantile(vector, c(0.75), na.rm = TRUE)
        iqr <- p75 - p25
        uplim <- p75 + 1.5*iqr
        lowlim <- p25 - 1.5*iqr
        numOutliers <- sum(vector < lowlim, na.rm = TRUE) + sum(vector > uplim, na.rm = TRUE)
        return(numOutliers)    
      }
    }
    varNames <- names(data)
    classVar <- sapply(data, class)
    conVars <- varNames[classVar == 'numeric' | classVar == 'integer']
    len <- length(conVars)
    if(len > 0) {
      if(class(data)[1] == 'data.frame') {
        dataConVar <- data[conVars]
      } else if (class(data)[1] == 'data.table'){
        dataConVar <- data[, conVars, with = FALSE]
      }
      variable <- names(dataConVar)
      nonMissingValues <- sapply(dataConVar, function(var) sum(!is.na(var)))
      missingValues <- sapply(dataConVar, function(var) sum(is.na(var)))
      missingPercentage <- sapply(missingValues, function(value) round(value/nrow(data)*100, digits = 2))
      mean <- sapply(dataConVar, function(var) mean(var, na.rm = TRUE))
      maximum <- sapply(dataConVar, function(var) max(var, na.rm = TRUE))
      minimum <- sapply(dataConVar, function(var) min(var, na.rm = TRUE))
      stdDeviation <-  sapply(dataConVar, function(var) sd(var, na.rm = TRUE))
      variance <- sapply(dataConVar, function(var) var(var, na.rm = TRUE))
      percentiles <- t(sapply(dataConVar, function(var) 
        quantile(var, c(.01, .05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, .99), na.rm = TRUE)))
      numOutliers <- sapply(dataConVar, numoutliers)
      conVarSummary <- data.frame(variable, 
                                  nonMissingValues, 
                                  missingValues,
                                  missingPercentage,
                                  minimum,
                                  mean, 
                                  maximum,
                                  stdDeviation,
                                  variance,
                                  percentiles,
                                  numOutliers,
                                  row.names = NULL)
      names(conVarSummary)[10:18] <- c('p01', 'p05', 'p10', 'p25', 'p50', 'p75', 'p90', 'p95', 'p99')
      return(conVarSummary)
    } else {
      stop('There are no variables of class integer or numeric in the data')
    }
  }
}