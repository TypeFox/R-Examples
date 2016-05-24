#' get mean rainfall bar plot of the input dataset or time series.
#' 
#' get mean rainfall bar plot of the input dataset or time series.
#' 
#' 
#' @param data A list containing different information, should be the result of reading netcdf file using
#' \code{\link{loadNcdf}}, or a time series, with first column the Date, second the value.
#' Time series can be an ENSEMBLE containning different members. Than the mean value will be given and the range will be given.
#' @param method A string showing the calculating method of the input time series. More information
#' please refer to the details.
#' @param cell A vector containing the locaton of the cell, e.g. c(2, 3), default is "mean", representing
#' the spatially averaged value. Check details for more information.
#' @param output A string showing the type of the output, if \code{output = 'ggplot'}, the returned 
#' data can be used in ggplot and \code{getPreciBar_comb()}; if \code{output = 'plot'}, the returned data is the plot containing all 
#' layers' information, and can be plot directly or used in grid.arrange; if not set, the data
#' will be returned.
#' @param name If \code{output = 'ggplot'}, name has to be assigned to your output, in order to differentiate
#' different outputs in the later multiplot using \code{getSpatialMap_comb}.
#' @param plotRange A boolean showing whether the range will be plotted.
#' @param member A number showing which member is selected to get, if the dataset has a "member" dimension. Default
#' is NULL, if no member assigned, and there is a "member" in dimensions, the mean value of the members will be
#' taken.
#' @param omitNA A boolean showing whether the missing value is omitted.
#' @param info A boolean showing whether the information of the map, e.g., max, mean ..., default is FALSE.
#' @param ... \code{title, x, y} showing the title and x and y axis of the plot. e.g. \code{title = 'aaa'}
#' @details
#' There are following methods to be selected, 
#' "annual": annual rainfall of each year is plotted.  
#' "winter", "spring", "autumn", "summer": seasonal rainfall of each year is plotted.
#' Month(number 1 to 12): month rainfall of each year is plotted, e.g. march rainfall of each year.
#' "meanMonthly": the mean monthly rainfall of each month over the whole period.
#' 
#' #Since "winter" is a crossing year, 12, 1, 2, 12 is in former year, and 1, 2 are in latter year.
#' #so winter belongs to the latter year.
#' 
#' 
#' \code{cell} representing the location of the cell, NOTE: this location means the index of the cell,
#' IT IS NOT THE LONGITUDE AND LATITUDE. e.g., \code{cell = c(2, 3)}, the program will take the 2nd longitude
#' and 3rd latitude, by the increasing order. Longitude comes first.
#' 
#' 
#' It is a generic function, if in your case you need to debug, please see \code{?debug()} 
#' for how to debug S4 method.
#' 
#' @examples
#' #gridData provided by package is the result of \code{loadNcdf()}
#' data(tgridData)
#' b1 <- getPreciBar(tgridData, method = 'annual')
#' b2 <- getPreciBar(tgridData, method = 'meanMonthly')
#' 
#' data(testdl)
#' TS  <- testdl[[1]]
#' a <- getPreciBar(TS, method = 'spring')
#' # if info = T, the information will be given at the bottom.
#' a <- getPreciBar(TS, method = 'spring', info = TRUE)
#' 
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @references 
#' 
#' 
#' \itemize{
#' \item Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software,
#' 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.
#' \item H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York, 2009.
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' }
#'
#' 
#' @return The calculated mean value of the input time series and the plot of the result.
#' @export
setGeneric('getPreciBar', function(data, method, cell = 'mean', output = 'data', name = NULL, 
                                   plotRange = TRUE, member = NULL, omitNA = TRUE, info = FALSE,
                                   ...) {
  standardGeneric('getPreciBar')
})

#' @describeIn getPreciBar
setMethod('getPreciBar', signature('list'), 
          function(data, method, cell, output, name, plotRange, member, omitNA, info, ...) {
            TS <- getPreciBar.list(data, cell, member)
            # for hyfo file, in order to process the data, year and month index need to be provided.
            startTime <- as.POSIXlt(data$Dates$start, tz = 'GMT')
            yearIndex <- startTime$year + 1900
            monthIndex <- startTime$mon + 1
            
            result <- getPreciBar.plot(TS, method, output, name, plotRange, omitNA, info, yearIndex,
                                       monthIndex, ...)
            return(result)
})

#' @describeIn getPreciBar
setMethod('getPreciBar', signature('data.frame'), 
          function(data, method, cell, output, name, plotRange, member, omitNA, info, ...) {
            Date <- as.POSIXlt(TS[, 1])
            yearIndex <- Date$year + 1900
            monthIndex <- Date$mon + 1
            TS <- getPreciBar.TS(data)
            result <- getPreciBar.plot(TS, method, output, name, plotRange, omitNA, info, 
                                       yearIndex, monthIndex, ...)
            return(result)
})


getPreciBar.list <- function(dataset, cell, member) {
  #check input dataset
  checkHyfo(dataset)
  
  data <- dataset$Data
  
  # Dimension needs to be arranged. Make sure first and second dimension is lat and lon.
  data <- adjustDim(data, ref = c('lon', 'lat', 'time'))
  
  # Because in the following part, only 3 dimensions are allowed, so data has to be processed.
  if (is.null(member) & any(attributes(data)$dimensions == 'member')) {
    dimIndex3 <- which(attributes(data)$dimensions != 'member')
    data <- apply(data, MARGIN = dimIndex3, FUN = mean, na.rm = TRUE)
  } else if (!is.null(member) & any(attributes(data)$dimensions == 'member')) {
    dimIndex3 <- which(attributes(data)$dimensions == 'member')
    data <- chooseDim(data, dimIndex3, member, drop = TRUE)
  } else if (!is.null(member) & !any(attributes(data)$dimensions == 'member')){
    stop('There is no member part in the dataset, but you choose one, check the input
         dataset or change your arguments.')
  }
  
  if (identical(cell, 'mean')) {
    TS <- apply(data, MARGIN = 3, FUN = mean, na.rm = TRUE) 
  } else {
    TS <- data[cell[1], cell[2], ]
  }
  
  return(TS)
}


#' @importFrom reshape2 melt
getPreciBar.TS <- function(TS) {
  
#  Date <- as.POSIXlt(TS[, 1])
#  yearIndex <- Date$year + 1900
#  monthIndex <- Date$mon + 1
  n <- ncol(TS) - 1
  
  if ( n == 1) {
    TS <- TS[, 2]
  } else {
    
    TS <- TS[, -1]
    # month index should be repeat, but years cannot.
#    yearIndex <- sapply(1:n, function(x) yearIndex + x - 1)
#    dim(yearIndex) <- c(n * nrow(yearIndex), 1)
    
#    monthIndex <- rep(monthIndex, n)
    TS <- melt(TS)[, 2]
    
  }
  return(TS)
}


#' @importFrom stats median
#' @importFrom reshape2 melt
#' @import ggplot2
getPreciBar.plot <- function(TS, method, output, name, plotRange, omitNA, info, 
                             yearIndex = NULL, monthIndex = NULL, ...) {
  
  
  if (method == 'meanMonthly') {
    
    monthlyPreci <- tapply(TS, INDEX = list(yearIndex, monthIndex), FUN = sum, na.rm = omitNA)
    meanMonthlyPreci <- apply(monthlyPreci, MARGIN = 2, FUN = mean, na.rm = TRUE)
    
    
    title <- 'Mean Monthly Precipitation'
    xlab <- 'Month'
    
    plotPreci <- data.frame(Index = month.abb[as.numeric(colnames(monthlyPreci))], 
                            Preci = meanMonthlyPreci)
    
    # Here factor has to be reassigned, to keep the original order, or it will be reordered.
    plotPreci$Index <- factor(plotPreci$Index, levels = plotPreci$Index, ordered = TRUE)
    
    if (plotRange) {
      maxValue <- apply(monthlyPreci, MARGIN = 2, FUN = max, na.rm = TRUE)
      minValue <- apply(monthlyPreci, MARGIN = 2, FUN = min, na.rm = TRUE)
      
      plotPreci$maxValue <- maxValue
      plotPreci$minValue <- minValue
      
      ylim <- c(0,max(maxValue, na.rm = TRUE) * 1.1)
      
    } else {
      ylim <- c(0,max(meanMonthlyPreci, na.rm = TRUE) * 1.1)
    }
    
    
  } else if (method == 'annual') {  
    
    if (length(unique(monthIndex)) < 12) {
      warning ('There are less than 12 months in a year, the results may be inaccurate.')
    }
    
    annualPreci <- tapply(TS, INDEX = yearIndex, FUN = sum, na.rm = TRUE)
    title <- 'Annual Precipitation'
    xlab <- 'Year'
    plotName <- names(annualPreci)
    
    plotPreci <- data.frame(Index = names(annualPreci), Preci = annualPreci)
    plotPreci$Index <- factor(plotPreci$Index, levels = plotPreci$Index, ordered = TRUE)
    
    ylim <- c(0, max(annualPreci, na.rm = TRUE) * 1.1)
    
  } else if (is.numeric(method)) {
    month <- method
    monExisting <- length(which(unique(monthIndex) == month))
    if (monExisting == 0) stop("Your input month doesn't exist in the dataset.")
    
    monthlyPreci <- getMeanPreci(TS, method = month, yearIndex = yearIndex,
                                 monthIndex = monthIndex, fullResults = TRUE, omitNA = omitNA)
    # If monthlyPreci length is 1, names need to be added.
    if (length(monthlyPreci) == 1) names(monthlyPreci) <- unique(yearIndex)
    plotPreci <- data.frame(Index = names(monthlyPreci), Preci = monthlyPreci)
    plotPreci$Index <- factor(plotPreci$Index, levels = plotPreci$Index, ordered = TRUE)
    
    title <- paste(month.abb[month], 'Precipitation over Whole Period', sep = ' ')
    xlab <- 'Year'
    ylim <- c(0, max(monthlyPreci, na.rm = TRUE) * 1.1)
    
  } else if (method == 'spring') {   
    
    wm <- match(c(3, 4, 5), unique(monthIndex))
    if (length(which(!is.na(wm))) < 3) {
      stop('Spring has less than 3 months, check data and try to calculate every month
           seperately or choose another season.')
    }
    
    seasonalPreci <- getMeanPreci(TS, method = 'spring', yearIndex = yearIndex,
                                  monthIndex = monthIndex, fullResults = TRUE, omitNA = omitNA)
    plotPreci <- data.frame(Index = names(seasonalPreci), Preci = seasonalPreci)
    plotPreci$Index <- factor(plotPreci$Index, levels = plotPreci$Index, ordered = TRUE)
    
    title <- paste('Spring', 'Precipitation over Whole Period', sep = ' ')
    xlab <- 'Year'
    ylim <- c(0, max(seasonalPreci, na.rm = TRUE) * 1.1)
    
    
    } else if (method == 'summer') {
      
      wm <- match(c(6, 7, 8), unique(monthIndex))
      if (length(which(!is.na(wm))) < 3) {
        stop('Summer has less than 3 months, check data and try to calculate every month
             seperately or choose another season.')
      }
      
      seasonalPreci <- getMeanPreci(TS, method = 'summer', yearIndex = yearIndex,
                                    monthIndex = monthIndex, fullResults = TRUE, omitNA = omitNA)
      plotPreci <- data.frame(Index = names(seasonalPreci), Preci = seasonalPreci)
      plotPreci$Index <- factor(plotPreci$Index, levels = plotPreci$Index, ordered = TRUE)
      
      title <- paste('Summer', 'Precipitation over Whole Period', sep = ' ')
      xlab <- 'Year'
      ylim <- c(0, max(seasonalPreci, na.rm = TRUE) * 1.1)
      
      
      } else if (method == 'autumn') {
        wm <- match(c(9, 10, 11), unique(monthIndex))
        if (length(which(!is.na(wm))) < 3) {
          stop('Autumn has less than 3 months, check data and try to calculate every month
               seperately or choose another season.')
        }
        
        seasonalPreci <- getMeanPreci(TS, method = 'autumn', yearIndex = yearIndex,
                                      monthIndex = monthIndex, fullResults = TRUE, omitNA = omitNA)
        plotPreci <- data.frame(Index = names(seasonalPreci), Preci = seasonalPreci)
        plotPreci$Index <- factor(plotPreci$Index, levels = plotPreci$Index, ordered = TRUE)
        
        title <- paste('Autumn', 'Precipitation over Whole Period', sep = ' ')
        xlab <- 'Year'
        ylim <- c(0, max(seasonalPreci, na.rm = TRUE) * 1.1)
        
        } else if (method == 'winter') {
          wm <- match(c(12, 1, 2), unique(monthIndex))
          if (length(which(!is.na(wm))) < 3) {
            stop('Winter has less than 3 months, check data and try to calculate every month
                 seperately or choose another season.')
          }
          
          seasonalPreci <- getMeanPreci(TS, method = 'winter', yearIndex = yearIndex,
                                        monthIndex = monthIndex, fullResults = TRUE, omitNA = omitNA)
          plotPreci <- data.frame(Index = names(seasonalPreci), Preci = seasonalPreci)
          plotPreci$Index <- factor(plotPreci$Index, levels = plotPreci$Index, ordered = TRUE)
          
          title <- paste('Winter', 'Precipitation over Whole Period', sep = ' ')
          xlab <- 'Year'
          ylim <- c(0, max(seasonalPreci, na.rm = TRUE) * 1.1)
          
          } else {
            stop(paste('No method called "', method, '", check help for information'))
          }
  
  
  xlim <- c(0, length(rownames(plotPreci))) 
  
  if (info == TRUE) {
    meanValue <- round(mean(plotPreci$Preci, na.rm = TRUE), 2)
    medianValue <- round(median(plotPreci$Preci,na.rm = TRUE), 2)
    plotMean <- paste('Mean', ' = ', meanValue)
    plotMedian <- paste('Median', ' = ', medianValue)
    
    plotMax <- round(max(plotPreci$Preci, na.rm = TRUE), 2)
    plotMin <- round(min(plotPreci$Preci, na.rm = TRUE), 2)
    word <- paste('\n\n', paste(' Max', '=', plotMax), ',', paste('Min', '=', plotMin), ',',
                  plotMean, ',', plotMedian)
  } else word <- NULL
  
  
  xlab <- paste(xlab, word)
  
  theme_set(theme_bw())
  
  mainLayer <- with(plotPreci, {
    ggplot(plotPreci) +
      geom_bar(aes(x = Index, y = Preci), stat = 'identity', colour = 'black', fill = 'cyan2', width = rel(.4)) +
      xlab(xlab) +
      ylab('Precipitation (mm)') +
      ggtitle(title) +
      labs(empty = NULL, ...) +#in order to pass "...", arguments shouldn't be empty.
      theme(plot.title = element_text(size = rel(1.6), face = 'bold'),
            axis.title.x = element_text(size = rel(1.6)),
            axis.title.y = element_text(size = rel(1.6)),
            axis.text.x = element_text(angle = 90, hjust = 1, size = rel(1.9)),
            axis.text.y = element_text(size = rel(1.9)))
    #    geom_text(x = min(xlim) + 0.95 * (max(xlim) - min(xlim)), y = min(ylim) + 0.15 * (max(ylim) - min(ylim)),
    #              label = word)+
    #     geom_hline(yintercept = meanValue) +
    #     geom_text(x = min(xlim) + 0.3 * (max(xlim) - min(xlim)), y = meanValue + 3, vjust = 0, label = 'mean') +
    #     geom_hline(yintercept = medianValue, colour = 'red') +
    #     geom_text(x = min(xlim) + 0.6 * (max(xlim) - min(xlim)), y = medianValue + 3, vjust = 0,
    #               label = 'median', colour = 'red')
  })
  
  
  if (plotRange) {
    if (is.null(plotPreci$maxValue)) {
      message('There is no plotRange for this method')
      print(mainLayer)
    } else {
      rangeLayer <- with(plotPreci, {
        geom_errorbar(aes(x = Index, ymax = maxValue, ymin = minValue), width = rel(0.3))
      })       
      print(mainLayer + rangeLayer)
    }
    
  } else {
    print(mainLayer)
  } 
  
  if (output == 'plot') {
    return(mainLayer)
  } else if (output == 'ggplot') {
    if (is.null(name)) stop('"name" argument not found, 
                            If you choose "ggplot" as output, please assign a name.')
    plotPreci$Name <- rep(name, dim(plotPreci)[1])
    return(plotPreci)
  } else {
    return(plotPreci)
  }
}







#' Combine bars together
#' @param ... different barplots generated by \code{getPreciBar(, output = 'ggplot')}, refer to details.
#' @details
#' ..., representing different ouput generated by \code{getPreciBar(, output = 'ggplot')}, they 
#' have to be of the same type, e.g., 
#' 1. Jan precipitation of different years, Feb precipitation of different years, and... 
#' They are both monthly precipitation, and they share x axis.
#' 
#' 2. Mean monthly precipitation of different dataset. e.g., long term mean monthly precipitation
#' and short term mean monthly precipitation. They are both mean monthly precipitation.
#' 
#' @param nrow A number showing the number of rows.
#' @param list If input is a list containing different ggplot data, use l\code{list = inputlist}.
#' NOTE: yOU HAVE TO PUT A \code{list = }, before your list.
#' @param x A string of x axis name.
#' @param y A string of y axis name.
#' @param title A string of the title.
#' @param output A boolean, if chosen TRUE, the output will be given.
#' @return A combined barplot.
#' @examples
#' 
#' data(tgridData)# the result of \code{\link{loadNcdf}}
#' #output type of getPreciBar() has to be 'ggplot'.
#' b1 <- getPreciBar(tgridData, method = 2, output = 'ggplot', name = 'b1')
#' b2 <- getPreciBar(tgridData, method = 3, output = 'ggplot', name = 'b2')
#' 
#' getPreciBar_comb(b1, b2)
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @export
#' @import ggplot2
#' @references 
#' 
#' \itemize{
#' \item H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York, 2009.
#' }
#' 
getPreciBar_comb <- function(..., list = NULL, nrow = 1, x = '', y = '', title = '', output = FALSE) {
  if (!is.null(list)) {
    data_ggplot <- do.call('rbind', list)
  } else {
    
    bars <- list(...)
    checkBind(bars, 'rbind')
    data_ggplot <- do.call('rbind', bars)
  }
  
  if (!class(data_ggplot) == 'data.frame') {
    warning('Your input is probably a list, but you forget to add "list = " before it.
            Try again, or check help for more information.')
  } else if (is.null(data_ggplot$Name)) {
    stop('No "Name" column in the input data, check the arguments in getPreciBar(), if 
         output = "ggplot" is assigned, more info please check ?getPreciBar.')
  }
  
  data_ggplot$Name <- factor(data_ggplot$Name, levels = unique(data_ggplot$Name), ordered = TRUE)
  
  theme_set(theme_bw())
  
  mainLayer <- with(data_ggplot, {
    ggplot(data_ggplot) +
      geom_bar(aes(x = Index, y = Preci),fill = 'cyan2', stat = 'identity', 
               colour = 'black', width = rel(.4)) +
      facet_wrap( ~ Name, nrow = nrow) +
      theme(plot.title = element_text(size = rel(1.6), face = 'bold'),
            axis.title.x = element_text(size = rel(1.6)),
            axis.title.y = element_text(size = rel(1.6)),
            axis.text.x = element_text(angle = 90, hjust = 1, size = rel(1.9)),
            axis.text.y = element_text(size = rel(1.9))) +
      labs(x = x, y = y, title = title)
  })
  
  if (!any(is.na(match(c('minValue', 'maxValue'), colnames(data_ggplot))))) {
    rangeLayer <- with(data_ggplot, {
      geom_errorbar(aes(x = Index, ymax = maxValue, ymin = minValue), width = rel(0.3))
    })       
    mainLayer <- mainLayer + rangeLayer
  }
  
  
  suppressWarnings(print(mainLayer))
  
  if (output == TRUE) return(data_ggplot)
  }

