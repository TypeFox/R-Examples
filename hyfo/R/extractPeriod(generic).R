#' Extract period from list or dataframe.
#' 
#' Extract common period or certain period from a list of different dataframes of time series, or from a 
#' dataframe. 
#' NOTE: all the dates in the datalist should follow the format in ?as.Date{base}.
#' @param data A list of different dataframes of time series, or a dataframe with first column Date, the rest columns value.
#' @param startDate A Date showing the start of the extract period, default as NULL, check details.
#' @param endDate A Date showing the end of the extract period, default as NULL, check details.
#' @param commonPeriod A boolean showing whether the common period is extracted. If chosen, startDate and endDate
#' should be NULL.
#' @param year extract certain year in the entire time series. if you want to extract year 2000, set \code{year = 2000}
#' @param month extract certain months in a year. e.g. if you want to extract Jan, Feb of each year, 
#' set \code{month = c(1, 2)}.
#' @details 
#' \strong{startDate and endDate}
#' 
#' If startDate and endDate are assigned, then certain period between startDate and endDate will be returned, 
#' for both datalist input and dataframe input.
#' 
#' If startDate and endDate are NOT assigned, then,
#' 
#'    if input is a datalist, the startDate and endDate of the common period of different datalists will be assigned
#'    to the startDate and endDate.
#' 
#'    if input is a dataframe, the startDate and endDate of the input dataframe will be assigned to the startDate
#'    and endDate . Since different value columns share a common Date column in a dataframe input. 
#' 
#' \strong{year and month}
#' 
#' For year crossing month input, hyfo will take from the year before. E.g. if \code{month = c(10, 11, 12, 1)},
#' and \code{year = 1999}, hyfo will take month 10, 11 and 12 from year 1998, and month 1 from 1999.You DO NOT 
#' have to set \code{year = 1998 : 1999}.
#' 
#' Well, if you set \code{year = 1998 : 1999}, hyfo will take month 10, 11 and 12 from year 1997, and month 1 from 1998,
#' then, take month 10, 11 and 12 from year 1998, month 1 from 1999. So you only have to care about the latter year.
#' 
#' It is a generic function, if in your case you need to debug, please see \code{?debug()} 
#' for how to debug S4 method.
#' 
#' @return A list or a dataframe with all the time series inside containing the same period.
#' @examples
#' # Generate timeseries datalist. Each data frame consists of a Date and a value.
#' 
#' AAA <- data.frame(
#' # date column
#' Date = seq(as.Date('1990-10-28'),as.Date('1997-4-1'),1),
#'  # value column
#' AAA = sample(1:100,length(seq(as.Date('1990-10-28'),as.Date('1997-4-1'),1)), repl = TRUE))
#' 
#' BBB <- data.frame(
#' Date = seq(as.Date('1993-3-28'),as.Date('1999-1-1'),1), 
#' BBB = sample(1:100,length(seq(as.Date('1993-3-28'),as.Date('1999-1-1'),1)), repl = TRUE))
#'  
#' CCC <- data.frame(
#' Date = seq(as.Date('1988-2-2'),as.Date('1996-1-1'),1), 
#' CCC = sample(1:100,length(seq(as.Date('1988-2-2'),as.Date('1996-1-1'),1)), repl = TRUE)) 
#' 
#' list <- list(AAA, BBB, CCC)# dput() and dget() can be used to save and load list file.
#' 
#' list_com <- extractPeriod(list, commonPeriod = TRUE)
#' 
#' # list_com is the extracted datalist.
#' str(list_com)
#' 
#' # If startDate and endDate is provided, the record between them will be extracted.
#' # make sure startDate is later than any startDate in each dataframe and endDate is 
#' # earlier than any endDate in each dataframe.
#' 
#' data(testdl)
#' datalist_com1 <- extractPeriod(testdl, startDate = '1994-1-1', endDate = '1995-10-1')
#' 
#' 
#' dataframe <- list2Dataframe(datalist_com1)
#' # now we have a dataframe to extract certain months and years.
#' dataframe_new <- extractPeriod(dataframe, month = c(1,2,3))
#' dataframe_new <- extractPeriod(dataframe, month = c(12,1,2), year = 1995)
#' 
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @importFrom zoo as.Date
#' @references 
#' 
#' \itemize{
#' \item Achim Zeileis and Gabor Grothendieck (2005). zoo: S3 Infrastructure for Regular and Irregular Time
#' Series. Journal of Statistical Software, 14(6), 1-27. URL http://www.jstatsoft.org/v14/i06/
#' }
#'
#' @export
setGeneric('extractPeriod', function(data, startDate = NULL, endDate = NULL, commonPeriod = FALSE, 
                                     year = NULL, month = NULL) {
  standardGeneric('extractPeriod')
})


#' @describeIn extractPeriod
#' @importFrom methods setMethod
setMethod('extractPeriod', signature('data.frame'),
          function(data, startDate, endDate, commonPeriod, year, month) {
            dataframe <- data
            dataset <- extractPeriod_dataframe(dataframe, startDate = startDate, endDate = endDate, year = year,
                                               month = month)
            return(dataset)
          
})


#' @describeIn extractPeriod
#' @importFrom methods setMethod
setMethod('extractPeriod', signature('list'),
          function(data, startDate, endDate, commonPeriod, year, month) {
            datalist <- data
            if (!is.null(startDate) & !is.null(endDate) & commonPeriod == FALSE) {
              dataset <- lapply(data, extractPeriod_dataframe, startDate = startDate, endDate = endDate, year = year,
                                month = month)
            } else if (is.null(startDate) & is.null(endDate) & commonPeriod == TRUE) {
              
              Dates <- lapply(datalist, extractPeriod_getDate) 
              Dates <- do.call('rbind', Dates)
              
              startDate <- as.Date(max(Dates[, 1]))
              endDate <- as.Date(min(Dates[, 2]))
              
              dataset <- lapply(datalist, extractPeriod_dataframe, startDate = startDate, endDate = endDate, year = year,
                                month = month)
              
            } else {
              stop('Enter startDate and endDate, set commonPeriod as False, or simply set commonPeriod as TRUE')
            }
            return(dataset)
          })




extractPeriod_dataframe <- function(dataframe, startDate, endDate, year = NULL, month = NULL) {
  # to check whether first column is a date format
  if (!grepl('-|/', dataframe[1, 1])) {
    stop('First column is not date or Wrong Date formate, check the format in ?as.Date{base} 
         and use as.Date to convert.')
  }    
  dataframe[, 1] <- as.Date(dataframe[, 1])
  
  if (is.null(startDate)) startDate <- dataframe[1, 1]
  if (is.null(endDate)) endDate <- tail(dataframe[, 1], 1)
  
  startIndex <- which(dataframe[, 1] == startDate)
  endIndex <- which(dataframe[, 1] == endDate)
  if (length(startIndex) == 0 | length(endIndex) == 0) {
    stop('startDate and endDate exceeds the date limits in dataframe. Check datalsit please.')
  }
  output <- dataframe[startIndex:endIndex, ]
  
  
  if (!is.null(year)) {
    Date <- as.POSIXlt(output[, 1])
    yea <- Date$year + 1900
    mon <- Date$mon + 1
    
    if (is.null(month) || !any(sort(month) != month)) {
      DateIndex <- which(yea %in% year)
      if (length(DateIndex) == 0) stop('No input years in the input ts, check your input.')
      
      output <- output[DateIndex, ]
      
      # if year crossing  than sort(month) != month, in this case we need to
      # take months from last year.
    } else {
      
      
      startIndex <- intersect(which(yea == year[1] - 1), which(mon == month[1]))[1]
      endIndex <- tail(intersect(which(yea == tail(year, 1)), which(mon == tail(month, 1))), 1)
      
      
      if (is.na(startIndex) || length(endIndex) == 0 || startIndex > endIndex) {
        stop('Cannot find input months and input years in the input time series.')
      }
      output <- output[startIndex:endIndex, ]
      
      if (any(diff(year) != 1)) {
        # if year is not continuous, like 1999, 2003, 2005, than we have to sift again.  
        Date <- as.POSIXlt(output[, 1])
        yea <- Date$year + 1900
        mon <- Date$mon + 1
        
        DateIndex <- unlist(sapply(year, function(x) {
          startIndex <- intersect(which(yea == x - 1), which(mon == month[1]))[1]
          endIndex <- tail(intersect(which(yea == x), which(mon == tail(month, 1))), 1)
          index <- startIndex:endIndex
          return(index)
        }))
        
        
        output <- output[DateIndex, ]
        
        # cannot directly return output here, because sometimes, month can be incontinuous,
        # we still need the next process to sift month.
      }
    }
    
  }
  
  
  if (!is.null(month)) {
    Date <- as.POSIXlt(output[, 1])
    mon <- Date$mon + 1
    
    # %in% can deal with multiple equalities
    DateIndex <- which(mon %in% month)
    
    if (length(DateIndex) == 0) stop('No input months in the input ts, check your input.')
    
    output <- output[DateIndex, ]
  }
  
  
  return(output)  
  }


#' @importFrom utils tail
#' @references 
#' 
#' \itemize{
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' }
#' 
#' 
extractPeriod_getDate <- function(dataset) {
  
  if (!grepl('-|/', dataset[1, 1])) {
    stop('First column is not date or Wrong Date formate, check the format in ?as.Date{base}, 
         and use as.Date to convert.')
  }
  start <- as.Date(dataset[1, 1])
  end <- as.Date(tail(dataset[, 1], 1))
  
  
  return(c(start, end))
  }



