#' Fill gaps in the rainfall time series.
#'
#' @param dataset A dataframe with first column the time, the rest columns are rainfall data of different gauges
#' @param corPeriod A string showing the period used in the correlation computing, 
#' e.g. daily, monthly, yearly.
#' @return The filled dataframe
#' @details
#' the gap filler follows the rules below:
#' 
#'  1. The correlation coefficient of every two columns (except time column) is calculated.
#' the correlation coefficient calculation can be based on 'daily', 'monthly', 'annual',
#' in each case, the daily data, the monthly mean daily data and annual mean daily data of 
#' each column will be taken in the correlation calculation.
#' 
#' Then the correlation matrix is got, then based on the matrix, for each column, 
#' the 1st, 2nd, 3rd,... correlated column will be got. So if there is missing value in the
#' column, it will get data from orderly 1st, 2nd, 3rd column.
#' 
#'  2. The  simple linear regress is calculated between every two columns. When generating the
#'  linear coefficient, the incept should be force to 0. i.e. y = a*x + b should be forec to 
#'  y = a*x.
#'  
#'  3. Gap filling. E.g., on a certain date, there is a missing value in column A, then the
#'  correlation order is column B, column C, column D, which means A should take values from
#'  B firstly, if B is also missing data, then C, then D.
#'  
#'  Assuming finally value from column C is taken. Then according to step 2, A = a*C, then the
#'  final value filled in column A is missing_in_A = a*value_in_C, a is the linear coeffcient.
#' 
#' @examples
#' b <- read.table(text = '        Date  AAA  BBB  CCC  DDD  EEE
#' 49 1999-12-15 24.8 21.4 25.6 35.0 17.4
#' 50 1999-12-16   NA  0.6  1.5  6.3  2.5
#' 51 1999-12-17   NA 16.3 20.3  NA 19.2
#' 52 1999-12-18   13  1.6 NA  6.3  0.0
#' 53 1999-12-19   10 36.4 12.5 26.8 24.9
#' 54 1999-12-20   NA  0.0  0.0  0.2  0.0
#' 55 1999-12-21  0.2  0.0  0.0  0.0  0.0
#' 56 1999-12-22  0.0  0.0  0.0  0.0  0.0')
#' 
#' b1 <- fillGap(b) # if corPeriod is missing, 'daily' is taken as default.
#' 
#' data(testdl)
#' a <- extractPeriod(testdl, commonPeriod = TRUE)
#' a1 <- list2Dataframe(a)
#' a2 <- fillGap(a1)
#' a3 <- fillGap(a1, corPeriod = 'monthly')
#' 
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @references
#' Gap fiiling method based on correlation and linear regression.
#' 
#' \itemize{
#' \item Hirsch, Robert M., et al. "Statistical analysis of hydrologic data." Handbook of hydrology. (1992): 17-1.
#' Salas, Jose D. "Analysis and modeling of hydrologic time series." Handbook of hydrology 19 (1993): 1-72.
#' 
#' }
#' 
#' 
#' @export
fillGap <- function(dataset, corPeriod = 'daily') {
  
  if (!grepl('-|/', dataset[1, 1])) {
    stop('First column is not date or Wrong Date formate, check the format in ?as.Date{base} 
          and use as.Date to convert.')
  }
  Date <- as.Date(dataset[, 1])
  data <- data.frame(dataset[, 2:dim(dataset)[2]])
  names <- colnames(data)
  
  corN <- fillGap_cor(data, corPeriod = corPeriod, Date = Date)
  cat('\nCorrelation Coefficient\n')
  print(corN)
  
  corOrder <- apply(corN, MARGIN = 1, FUN = function(x) order(-x))
  corOrder <- corOrder[2:dim(corOrder)[1], ]
  corOrderName <- t(apply(corOrder, MARGIN = 2, FUN = function(x) names[x]))
  
  cat('\nCorrelation Order\n')
  colnames(corOrderName) <- seq(1 : dim(corOrderName)[2])
  print(corOrderName)
  
  lmCoef <- fillGap_lmCoef(data, corOrder)
  cat('\nLinear Coefficients\n')
  rownames(lmCoef) <- seq(1 : dim(corOrderName)[2])
  print(t(lmCoef))
  
  output <- lapply(1:dim(data)[2], fillGap_column, data = data,
                   corOrder = corOrder, lmCoef = lmCoef)
  output <- data.frame(output)
  colnames(output) <- names
  
  output <- cbind(Date, output)
  
  return(output)
}


#' Get monthly rainfall
#' 
#' @param TS A rainfall time series.
#' @param year A list showing the year index of the time series.
#' @param mon A list showing the mon index of the time series.
#' @return the monthly rainfall matrix of the rainfall time series.
monthlyPreci <- function(TS, year, mon) {
  
  # monthly daily mean is used in order not to affected by missing values.
  monTS <- tapply(TS, INDEX = list(year, mon), FUN = mean, na.rm = TRUE)
  output <- t(monTS)
  dim(output) <- c(dim(monTS)[1] * dim(monTS)[2], 1)
  return(output)
}


fillGap_column <- function(i, data, corOrder, lmCoef) {
  TS <- data[, i] # extract target column
  l <- dim(data)[2] # length
  
  for (j in 1:l) {
    if (!any(is.na(TS))) break
    NAindex <- which(is.na(TS))
    TS[NAindex] <- round(lmCoef[j, i] * data[NAindex, corOrder[j, i]], 3)
    
    if (j == l) stop('Error: One time consists of all NA values')
  }
  
  return(TS)
}


#' @importFrom stats cor na.omit
#' @references 
#' 
#' \itemize{
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' }
#' 
#' 

fillGap_cor <- function(data, corPeriod = 'daily', Date) {
  
  names <- colnames(data)
  year <- format(Date, '%Y')
  
  if (corPeriod == 'monthly') {
    #based on monthly rainfall
    mon <- format(Date, '%m')
    monthlyPreci <- lapply(data, FUN = monthlyPreci, year = year, mon = mon)
    corData <- do.call('cbind', monthlyPreci)
  } else if (corPeriod == 'yearly') {
    year <- format(Date, '%Y')
    # yearly daily mean is used in order not to affected by missing values.
    annualPreci <- lapply(data, FUN = function(x) tapply(x, INDEX = year, FUN = mean, na.rm = TRUE))
    corData <- do.call('cbind', annualPreci)
  } else if (corPeriod == 'daily') {
    corData <- data
  } else {
    stop('Pleas choose among "daily", "monthly", "yearly".')
  }
  
  corData <- data.frame(na.omit(corData))
  colnames(corData) <- names
  
  corN <- cor(corData)
  
  return(corN)
  
} 

#' @importFrom utils combn
#' @importFrom stats coef lm
#' @references 
#' R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' 
fillGap_lmCoef <- function(data, corOrder) {
  l <- dim(data)[2]
  m <- diag(l)# m is the coeficients matrix
  m[lower.tri(m)] <- combn(data, 2, function(x) coef(lm(x[, 2] ~ x[, 1] + 0)))
  tm <- t(m)
  
  tm[lower.tri(tm)] <- combn(data, 2, function(x) coef(lm(x[, 1] ~ x[, 2] + 0)))
  
  m <- t(tm)
  
  lmCoef <- lapply(1 : l, function(x) m[x,corOrder[, x]])
  lmCoef <- do.call('rbind', lmCoef)
  rownames(lmCoef) <- colnames(data)
  
  return(t(lmCoef))
}

