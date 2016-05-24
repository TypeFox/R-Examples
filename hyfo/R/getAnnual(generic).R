#' Get annual rainfall of different rainfall time series
#' 
#' Get annual rainfall of different raninfall time series.
#' 
#' 
#' @param data A list containing different time series of different rainfall gauges. Or a dataframe with first column Date and the rest columns the value of different
#' gauging stations. Usually an output of \code{list2Dataframe}.
#' @param output A string showing the output output.
#' @param minRecords A number showing the minimum accept record number, e.g. for a normal 
#' year(365 days), if \code{minRecords = 360}, it means if a year has less than 360 records
#' of a year, it will be ignored in the mean annual value calculation. Only valid 
#' when \code{output = "mean"}, default is 355.
#' @param ... \code{title, x, y} showing the title and x and y axis of the plot. e.g. \code{title = 'aaa'}
#' @return The annual rainfall and the number of missing data of each year and each rainfall gauge, which 
#' will also be plotted. If output "mean" is seleted, the mean annual rainfall will be returned.
#' @details 
#' It is a generic function, if in your case you need to debug, please see \code{?debug()} 
#' for how to debug S4 method.
#' 
#' @examples
#' #datalist is provided by the package as a test.
#' data(testdl)
#' a <- getAnnual(testdl)
#' #set minRecords to control the calculation of annual rainfall.
#' b <- getAnnual(testdl, output = 'mean', minRecords = 350)
#' c <- getAnnual(testdl, output = 'mean', minRecords = 365)
#' 
#' a1 <- extractPeriod(testdl, comm = TRUE)
#' a2 <- list2Dataframe(a1)
#' getAnnual(a2)
#' 
#' a3 <- fillGap(a2)
#' getAnnual(a3)
#' 
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @export
#' @importFrom methods setGeneric
#' 
#' @references 
#' \itemize{
#' \item H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York, 2009.
#' \item Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software,
#' 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' }
#' 
#' 
setGeneric('getAnnual', function(data, output = 'series', minRecords = 355, 
                                 ...) {
  standardGeneric('getAnnual')
})

#' @describeIn getAnnual
#' @importFrom methods setMethod
setMethod('getAnnual', signature('data.frame'), 
          function(data, output, minRecords, ...) {
            result <- getAnnual.TS(data)
            getAnnual.plot(result, output, minRecords, ...)
            return(result)
})

#' @describeIn getAnnual
#' @importFrom methods setMethod
setMethod('getAnnual', signature('list'),
          function(data, output, minRecords, ...) {
            result <- getAnnual.list(data)
            getAnnual.plot(result, output, minRecords, ...)
            return(result)
          })

getAnnual.TS <- function(dataframe) {
  Date <- as.POSIXlt(dataframe[, 1])
  # Calculate how many gauging stations.
  stations <- colnames(dataframe)[2:ncol(dataframe)]
  
  data <- lapply(stations, function(x) {
    dataframe_new <- data.frame(Date, dataframe[, x])
    colnames(dataframe_new)[2] <- x
    getAnnual_dataframe(dataframe_new)
  })
  
  data <- do.call('rbind', data)
  #  After rbind, factor level has to be reassigned in order to be well plotted.
  data$Year <- factor(data$Year, levels = sort(unique(data$Year)), ordered = TRUE)
  rownames(data) <- NULL
  
  return(data)
}


getAnnual.list <- function(datalist) {
  data <- lapply(datalist, FUN = getAnnual_dataframe)
  data <- do.call('rbind', data)
  #  After rbind, factor level has to be reassigned in order to be well plotted.
  data$Year <- factor(data$Year, levels = sort(unique(data$Year)), ordered = TRUE)
  rownames(data) <- NULL
  return(data)
}

#' @import ggplot2 
#' @importFrom reshape2 melt
#' @importFrom stats aggregate
getAnnual.plot <- function(data, output, minRecords, ...) {
  theme_set(theme_bw())
  
  if (output == 'mean') {
    validData <- data[data$recordNum >= minRecords,]
    
    data <- aggregate(validData$AnnualPreci, list(validData$Name), mean)
    colnames(data) <- c('Name', 'AnnualPreci')
    
    mainLayer <- with(data, {
      ggplot(data)+
        geom_bar(aes(x = Name, y = AnnualPreci, fill = Name), stat = 'identity')+
        labs(empty = NULL, ...)#in order to pass "...", arguments shouldn't be empty.
      
    })
    
    print(mainLayer)
    
  } else {
    
    plotData <- with(data, {
      subset(data, select = c(Year, Name, NANum, AnnualPreci))
    })
    
    plotData <- melt(plotData, var.id = c('Year', 'Name'))
    
    
    mainLayer <- with(plotData, {
      ggplot(plotData) +
        geom_bar(aes(x = Year, y = value , fill = Name), 
                 stat = 'identity') +
        facet_grid(variable ~ Name, scale = 'free') +
        xlab('Year') +
        ylab(NULL) +
        labs(empty = NULL, ...) +#in order to pass "...", arguments shouldn't be empty.
        theme(plot.title = element_text(size = 20, face = 'bold', vjust = 1)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = rel(1.5)),
              axis.text.y = element_text(size = rel(1.5)))
      #      grid.arrange(mainLayer, ncol = 4)
      
    })
    
    
    print(mainLayer)
  }  
}








#' Get annual rainfall of the input time series.
#' 
#' @param dataset A dataframe containing one time series, e.g., rainfall from one gauging station.
#' the time should follow the format : "1990-1-1"
#' @return The annual rainfall of each year of the input station.
# @examples
# data(testdl)
# getAnnual_dataframe(testdl[[1]])
#' 
getAnnual_dataframe <- function(dataset) {
  
  if (!grepl('-|/', dataset[1, 1])) {
    stop ('First column is not date or Wrong Date formate, check the format in ?as.Date{base},
          and use as.Date to convert.')
  }
  Date <- as.Date(dataset[, 1])
  year <- format(Date, '%Y')
  yearUnique <- unique(year)
  #  yearUnique <- factor(yearUnique, levels = yearUnique, ordered = TRUE)
  calcuNum <- c(1:length(yearUnique))
  
  
  annualPreci <- tapply(dataset[, 2], INDEX = year, FUN = sum, na.rm = TRUE)
  recordNum <- tapply(dataset[, 2], INDEX = year, function(x) length(which(!is.na(x))))
  NANum <- tapply(dataset[, 2], INDEX = year, function(x) length(which(is.na(x))))
  
  
  name <- rep(colnames(dataset)[2], length(calcuNum))
  output <- data.frame(Year = as.numeric(yearUnique), Name = name, AnnualPreci = annualPreci,
                       recordNum, NANum)
  
  #output$Year <- factor(output$Year, levels = output$Year, ordered = TRUE)
  return(output)
}

