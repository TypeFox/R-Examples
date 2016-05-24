#' Get ensemble forecast from historical data.
#' 
#' getHisEnsem use historical data as the forecasting input time series.
#' 
#' @param TS A time series dataframe, with first column Date, and second column value.
#' @param example A vector containing two strings showing the start and end date, which represent the 
#' forecasting period. Check details for more information.
#'
#' the program will extract every possible period in TS you provided to generate the ensemble. Check details for 
#' more information.
#' @param interval A number representing the interval of each ensemble member. NOTE: "interval" takes
#' 365 as a year, and 30 as a month, regardless of leap year and months with 31 days. So if you want the interval 
#' to be 2 years, set \code{interval = 730}, which equals 2 * 365 ; if two months, set \code{interval = 60}; 
#' 2 days, \code{interval = 2}, for other numbers that cannot be divided by 365 or 30 without remainder, it will treat the 
#' number as days.By defualt interval is set to be 365, a year.
#' @param buffer A number showing how many days are used as buffer period for models. Check details for more
#' information.
#' 
#' @param plot A string showing whether the plot will be shown, e.g., 'norm' means normal plot (without any process), 
#' 'cum' means cummulative plot, default is 'norm'. For other words there will be no plot.
#' @param output A string showing which type of output you want. Default is "data", if "ggplot", the 
#' data that can be directly plotted by ggplot2 will be returned, which is easier for you to make series
#' plots afterwards. NOTE: If \code{output = 'ggplot'}, the missing value in the data will
#' be replaced by \code{mv}, if assigned, default mv is 0.
#' 
#' @param name If \code{output = 'ggplot'}, name has to be assigned to your output, in order to differentiate
#' different outputs in the later multiplot using \code{getEnsem_comb}.
#' 
#' @param mv A number showing representing the missing value. When calculating the cumulative value, 
#' missing value will be replaced by mv, default is 0.
#' @param ... \code{title, x, y} showing the title and x and y axis of the plot. e.g. \code{title = 'aaa'}
#' 
#' @details 
#' 
#' \code{example} E.g., if you have a time series from 2000 to 2010. Assuming you are in 2003,
#' you want to forecast the period from 2003-2-1 to 2003-4-1. Then for each year in your input
#' time series, every year from 1st Feb to 1st Apr will be extracted to generate the ensemble
#' forecasts. In this case your input example should be \code{example = c('2003-2-1', '2003-4-1')}
#' 
#' \code{interval} doesn't care about leap year and the months with 31 days, it will take 365 as a year, and 30 as a month.
#' e.g., if the interval is from 1999-2-1 to 1999-3-1, you should just set interval to 30, although the real interval is 28
#' days.
#' 
#' \code{example} and \code{interval} controls how the ensemble will be generated. e.g. if the time series is from 
#' 1990-1-1 to 2001-1-1.
#' 
#' if \code{example = c('1992-3-1', '1994-1-1')} and \code{interval = 1095}, note, 1095 = 365 * 3, so the program treat
#' this as 3 years.
#' 
#' Then you are supposed to get the ensemble consisting of following part:
#' 
#' 1. 1992-3-1 to 1994-1-1 first one is the example, and it's NOT start from 1990-3-1.
#' 2. 1995-3-1 to 1997-1-1 second one starts from 1993, because "interval" is 3 years.
#' 3. 1998-3-1 to 2000-1-1
#' 
#' because the last one "2000-3-1 to 2002-1-1", 2002 exceeds the original TS range, so it will not be included.
#' 
#' Sometimes, there are leap years and months with 31 days included in some ensemble part, in which case the length of the data will
#' be different, e.g., 1999-1-1 to 1999-3-1 is 1 day less than 2000-1-1 to 2000-3-1. In this situation,
#' the data will use example as a standard. If the example is 1999-1-1 to 1999-3-1, then the latter one
#' will be changed to 2001-1-1 to 2000-2-29, which keeps the start Date and change the end Date.
#' 
#' If the end date is so important that cannot be changed, try to solve this problem by resetting
#' the example period, to make the event included in the example.
#' 
#' Good set of example and interval can generate good ensemble.
#' 
#' \code{buffer}
#' Sometimes the model needs to run for a few days to warm up, before the forecast. E.g., if a forecast starts at
#' '1990-1-20', for some model like MIKE NAM model, the run needs to be started about 14 days. So the input timeseries
#' should start from '1990-1-6'.
#' 
#' Buffer is mainly used for the model hotstart. Sometimes the hot start file cannot contain all the parameters needed,
#' only some important parameters. In this case, the model needs to run for some time, to make other parameters ready
#' for the simulation.
#' 
#' 
#' \code{name}
#' Assuming you have two ggplot outputs, you want to plot them together. In this situation, you
#' need a name column to differentiate one ggplot output from the other. You can assigne this name
#' by the argument directly, name has to be assigned if \code{output = 'ggplot'} is selected,
#' @return A ensemble time series using historical data as forecast.
#' 
#' @examples
#' 
#' data(testdl)
#' 
#' a <- testdl[[1]]
#' 
#' # Choose example from "1994-2-4" to "1996-1-4"
#' b <- getHisEnsem(a, example = c('1994-2-4', '1996-1-4'))
#' 
#' # Default interval is one year, can be set to other values, check help for information.
#' 
#' # Take 7 months as interval
#' b <- getHisEnsem(a, example = c('1994-2-4', '1996-1-4'), interval = 210, plot = 'cum') 
#' # Take 30 days as buffer
#' b <- getHisEnsem(a, example = c('1994-2-4', '1996-1-4'), interval = 210, buffer = 30)
#' 
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' 
#' @importFrom reshape2 melt 
#' @importFrom grDevices rainbow
#' @import ggplot2
#' @references 
#' 
#' \itemize{
#' \item Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software,
#' 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.
#' \item H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York, 2009.
#' }
#' 
#' 
#' @export

getHisEnsem <- function (TS, example, interval = 365, buffer = 0, plot = 'norm', output = 'data', 
                         name = NULL, mv = 0, ...) {
  if (!grepl('-|/', TS[1, 1])) {
    stop('First column is not date or Wrong Date formate, check the format in ?as.Date{base} 
         and use as.Date to convert.')
  } else if (!grepl('-|/', example[1]) | !grepl('-|/', example[1])) {
    stop('Wrong date format in the example, check the format in ?as.Date{base} 
         and use as.Date to convert.')
  } else {
    
    
    
    TS[, 1] <- as.Date(TS[, 1])
    example <- as.Date(example ,tz = '')
    exL <- example[2] - example[1]
    # Test if example is in the range of the TS
    a <- which(TS[, 1] == example[1] | TS[, 1] == example[2])
    if (length(a) < 2) stop('Example is out of the time series, reset example.')
    
    
    
    if (interval %% 365 == 0) {
      d <- interval / 365
      
      # Get sequence of start and end date.
      
      startDate <- rev(seq(from = example[1], to = min(TS[, 1]), by = paste(-d, 'years')))
      endDate <- seq(from = example[2], to = max(TS[, 1]), by = paste(d, 'years'))

      n <- length(startDate) + length(endDate) - 1 # example is counted twice, should be subtracted.      
      
      # Generate full start date series.
      startDate <- seq(min(startDate), length = n, by = paste(d, 'years'))
      endDate <- startDate + exL
      
    } else if (interval %% 30) {
      d <- interval / 30
      
      # Get sequence of start and end date.
      
      startDate <- rev(seq(from = example[1], to = min(TS[, 1]), by = paste(-d, 'months')))
      endDate <- seq(from = example[2], to = max(TS[, 1]), by = paste(d, 'months'))
      
      n <- length(startDate) + length(endDate) - 1
      
      startDate <- seq(min(startDate), length = n, by = paste(d, 'months'))
      endDate <- startDate + exL
      
    } else {
      d <- interval
      
      # Get sequence of start and end date.
      
      startDate <- rev(seq(from = example[1], to = min(TS[, 1]), by = paste(-d, 'days')))
      endDate <- seq(from = example[2], to = max(TS[, 1]), by = paste(d, 'days'))
      
      n <- length(startDate) + length(endDate) - 1
      
      startDate <- seq(min(startDate), length = n, by = paste(d, 'days'))
      endDate <- startDate + exL
    }
    
    data <- mapply(FUN = function(x, y) extractPeriod_dataframe(dataframe = TS, startDate = x, endDate = y),
                   x = startDate, y = endDate)
    
    data <- lapply(1:n, function(x) data.frame(data[, x]))
    
    if (buffer > 0) {
      bufferStart <- example[1] - buffer
      bufferEnd <- example[1] - 1
      bufferTS <- extractPeriod_dataframe(TS, bufferStart, bufferEnd)
      
      data <- lapply(data, function(x) rbind(bufferTS, x))
      
    } else if (buffer < 0) {
      stop ('Buffer should be positive, or reset example.')
    }
    
    
    data_output <- list2Dataframe(data)
    colnames(data_output) <- c('Date', as.character(startDate))
    
    # Rearrange dataframe to make example the first column.
    ind <- match(c('Date', as.character(example[1])), colnames(data_output))
    # when use cbind, to ensure the output is also a dataframe, one inside cbind should be dataframe
    # Even output is alread a dataframe, but when ind is a single number, then output[ind] will
    # not be a dataframe, but an array.
    data_output <- cbind(data.frame(data_output[ind]), data_output[-ind])
    ex_date <- seq(from = example[1] - buffer, to = example[2], by = 1)
    data_output$Date <- ex_date
    colnames(data_output)[2] <- 'Observation'
    
    meanV <- apply(data_output[, 2:ncol(data_output)], MARGIN = 1, FUN = mean, na.rm = TRUE)
    
    data_output <- cbind(data.frame(Date = data_output[, 1]), Mean = meanV, 
                         data_output[, 2:ncol(data_output)])
    
    data_ggplot <- melt(data_output, id.var = 'Date')
    NAIndex <- is.na(data_ggplot$value)
    data_ggplot$nav <- rep(0, nrow(data_ggplot))
    data_ggplot$nav[NAIndex] <- 1
    
    if (plot == 'norm') {
      data_ggplot$value[NAIndex] <- mv
      
    } else if (plot == 'cum') {
      data_output[is.na(data_output)] <- mv
      cum <- cbind(data.frame(Date = data_output$Date), cumsum(data_output[2:ncol(data_output)]))
        
      data_ggplot <- melt(cum, id.var = 'Date')
    } else {
      stop('plot can only be "norm" or "cum", do not assign other words')
    }
    
    #generate different colors 
    colors = c('brown1', 'dodgerblue3', rainbow(n = length(unique(data_ggplot$variable)) - 2,
                                               start = 0.1))
    
    theme_set(theme_bw())
    mainLayer <- with(data_ggplot, {
      ggplot(data = data_ggplot) +
        aes(x = Date, y = value, color = variable, group = variable) +
        geom_line(size = 0.5) +
        geom_line(data = data_ggplot[data_ggplot$variable == 'Observation', ], size = 1.6) +
        geom_line(data = data_ggplot[data_ggplot$variable == 'Mean', ], size = 1.6) +
        geom_point(data = data_ggplot[NAIndex, ], size = 3, shape = 4, color = 'black') +
        scale_color_manual(values = colors) +
        labs(empty = NULL, ...) +#in order to pass "...", arguments shouldn't be empty.
        theme(axis.text.x = element_text(size = rel(1.8)),
              axis.text.y = element_text(size = rel(1.8)),
              axis.title.x = element_text(size = rel(1.8)),
              axis.title.y = element_text(size = rel(1.8)))
    })
    print(mainLayer)
    
    if (output == 'ggplot') {
      if (is.null(name)) stop('"name" argument not found, 
                            If you choose "ggplot" as output, please assign a name.')
      data_ggplot$name <- rep(name, nrow(data_ggplot))       
      data_ggplot$nav <- rep(0, nrow(data_ggplot))
      data_ggplot$nav[NAIndex] <- 1

      return(data_ggplot)
    } else {
      return(data_output)
    }
  }
}






#' Extract time series from forecasting data.
#' 
#' getFrcEnsem extract timeseries from forecasting data, if forecasting data has a member session
#' an ensemble time sereis will be returned, if forecasting data doesn't have a member session, a singe time
#' series will be returned.
#' 
#' @param dataset A list containing different information, should be the result of \code{\link{loadNcdf}}
#' @param cell A vector containing the locaton of the cell, e.g. c(2, 3), default is "mean", representing
#' the spatially averaged value. Check details for more information.
#' @param plot A string showing whether the plot will be shown, e.g., 'norm' means normal plot (without any process), 
#' 'cum' means cummulative plot, default is 'norm'. For other words there will be no plot.
#' @param output A string showing which type of output you want. Default is "data", if "ggplot", the 
#' data that can be directly plotted by ggplot2 will be returned, which is easier for you to make series
#' plots afterwards. NOTE: If \code{output = 'ggplot'}, the missing value in the data will
#' be replaced by \code{mv}, if assigned, default mv is 0.
#' @param name If \code{output = 'ggplot'}, name has to be assigned to your output, in order to differentiate
#' different outputs in the later multiplot using \code{getEnsem_comb}.
#' @param mv A number showing representing the missing value. When calculating the cumulative value, 
#' missing value will be replaced by mv, default is 0.
#' @param coord A coordinate of longitude and latitude. e.g. corrd = c(lon, lat). If coord is assigned,
#' cell argument will no longer be used.
#' @param ... \code{title, x, y} showing the title and x and y axis of the plot. e.g. \code{title = 'aaa'}
#' 
#' @details 
#' 
#' \code{cell} representing the location of the cell, NOTE: this location means the index of the cell,
#' IT IS NOT THE LONGITUDE AND LATITUDE. e.g., \code{cell = c(2, 3)}, the program will take the 2nd longitude
#' and 3rd latitude, by the increasing order. Longitude comes first.
#' 
#' \code{name}
#' Assuming you have two ggplot outputs, you want to plot them together. In this situation, you
#' need a name column to differentiate one ggplot output from the other. You can assigne this name
#' by the argument directly, If name is not assigned and \code{output = 'ggplot'} is selected, then
#' the system time will be selected as name column.
#' 
#' @examples 
#' 
#' filePath <- system.file("extdata", "tnc.nc", package = "hyfo")

#' # Then if you don't know the variable name, you can use \code{getNcdfVar} to get variable name
#' varname <- getNcdfVar(filePath)
#' nc <- loadNcdf(filePath, varname)
#' a <- getFrcEnsem(nc)
#' 
#' # If there is no member session in the dataset, a single time sereis will be extracted.
#' a1 <- getFrcEnsem(tgridData)
#' 
#' 
#' # The default output is spatially averaged, if there are more than one cells in the dataset, 
#' # the mean value of the cells will be calculated. While if you are interested in special cell, 
#' # you can assign the cell value. You can also directly use longitude and latitude to extract 
#' # time series.
#' 
#' getSpatialMap(nc, 'mean')
#' a <- getFrcEnsem(nc, cell = c(6,2))
#' 
#' # From the map, cell = c(6, 2) means lon = -1.4, lat = 43.2, so you can use corrd to locate
#' # your research area and extract time series.
#' b <- getFrcEnsem(nc, coord = c(-1.4, 43.2))
#' 
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @return A ensemble time series extracted from forecating data.
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @references 
#' 
#' \itemize{
#' \item H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York, 2009.
#' \item Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software,
#' 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.
#' \item Santander Meteorology Group (2015). downscaleR: Climate data manipulation and
#' statistical downscaling. R package version 0.6-0.
#' https://github.com/SantanderMetGroup/downscaleR/wiki
#' }
#' 
#' 
#' @export
getFrcEnsem <- function(dataset, cell = 'mean', plot = 'norm', output = 'data', name = NULL,
                        mv = 0, coord = NULL, ...) {
  # cell should be a vector showing the location, or mean representing the loacation averaged.
  
  checkHyfo(dataset)
  
  Date <- as.Date(dataset$Dates$start)
  data <- dataset$Data
  
  # Dimension needs to be arranged. Make sure first and second dimension is lat and lon.
  data <- adjustDim(data, ref = c('lon', 'lat', 'time'))
  
  if (!is.null(coord)) {
    cell <- coord2cell(coord, dataset$xyCoords$x, dataset$xyCoords$y)
  } 
  
  
  if (!any(attributes(data)$dimensions == 'member')){
    message('There is no member part in the dataset, there will be only one column of value
            returned.')
    
    if (length(cell) == 2) {
      data_ensem <- data[cell[1], cell[2], ]
      
    } else if (cell == 'mean') {
      data_ensem <- apply(data, MARGIN = 3, FUN = mean, na.rm = TRUE)
      #    colnames <- 1:ncol(data_ensem)
      
    } else {
      stop('Wrong cell input, check help for information.')
    }
    
  } else {
    
    if (length(cell) == 2) {
      data_ensem <- data[cell[1], cell[2], , ]
      meanV <- apply(data_ensem, MARGIN = 1, FUN = mean, na.rm = TRUE)
      data_ensem <- data.frame('Mean' = meanV, data_ensem) 
      
    } else if (cell == 'mean') {
      data_ensem <- apply(data, MARGIN = c(3, 4), FUN = mean, na.rm = TRUE)
  #    colnames <- 1:ncol(data_ensem)
      meanV <- apply(data_ensem, MARGIN = 1, FUN = mean, na.rm = TRUE)
      data_ensem <- data.frame('Mean' = meanV, data_ensem)
      
    } else {
      stop('Wrong cell input, check help for information.')
    }
  }

  
  data_output <- data.frame(Date, data_ensem)
  data_ggplot <- melt(data_output, id.var = 'Date')
  NAIndex <- is.na(data_ggplot$value)
  
  
  if (plot == 'norm') {
    data_ggplot$value[NAIndex] <- mv
  } else if (plot == 'cum') {
    data_output[is.na(data_output)] <- mv
    cum <- cbind(data.frame(Date = data_output$Date), cumsum(data_output[2:ncol(data_output)]))
    
    data_ggplot <- melt(cum, id.var = 'Date')
    
  }
  
  colors = c('brown1', rainbow(n = length(unique(data_ggplot$variable)) - 1,
                                              start = 0.1))
  
  theme_set(theme_bw())
  mainLayer <- with(data_ggplot, {
    ggplot(data = data_ggplot) +
      aes(x = Date, y = value, color = variable) +
      geom_line(size = 0.5) +
      geom_line(data = data_ggplot[data_ggplot$variable == 'Mean', ], size = 1.6, color = 'red') +
      geom_point(data = data_ggplot[NAIndex, ], size = 2, shape = 4, color = 'black') +
      scale_color_manual(values = colors) +
      theme(axis.text.x = element_text(size = rel(1.8)),
            axis.text.y = element_text(size = rel(1.8)),
            axis.title.x = element_text(size = rel(1.8)),
            axis.title.y = element_text(size = rel(1.8))) +
      labs(empty = NULL, ...)#in order to pass "...", arguments shouldn't be empty.
    
  })
  print(mainLayer)
  
  if (output == 'ggplot') {
    if (is.null(name)) stop('"name" argument not found, 
                            If you choose "ggplot" as output, please assign a name.')
    
    data_ggplot$name <- rep(name, nrow(data_ggplot))     
    data_ggplot$nav <- rep(0, nrow(data_ggplot))
    data_ggplot$nav[NAIndex] <- 1
    return(data_ggplot)
  } else {
    return(data_output)
  }
}



#' Combine ensembles together
#' @param ... different ensembles generated by \code{getHisEnsem(, output = 'ggplot')} 
#' or \code{getFrcEnsem(, output = 'ggplot')}, see details.
#' @param nrow A number showing the number of rows.
#' @param list If input is a list containing different ggplot data, use \code{list = inputlist}.
#' @param legend A boolean representing whether you want the legend. Sometimes when you combine
#' plots, there will be a lot of legends, if you don't like it, you can turn it off by setting
#' \code{legend = FALSE}, default is TRUE.
#' @param x A string of x axis name.
#' @param y A string of y axis name.
#' @param title A string of the title.
#' @param output A boolean, if chosen TRUE, the output will be given.
#' @return A combined ensemble plot.
#' @examples 
#' 
#' data(testdl)
#' 
#' a <- testdl[[1]]
#' 
#' # Choose example from "1994-2-4" to "1996-1-4"
#' 
#' 
#' b1<- getHisEnsem(a, example = c('1995-2-4', '1996-1-4'), plot = 'cum', output = 'ggplot',
#'                  name = 1)
#'                   
#' b2 <- getHisEnsem(a, example = c('1995-4-4', '1996-3-4'), plot = 'cum', output = 'ggplot',
#'                  name = 2)
#' 
#' getEnsem_comb(b1, b2)
#' getEnsem_comb(list = list(b1, b2), nrow = 2)
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' 
#' @export
#' @import ggplot2
#' @references 
#' 
#' \itemize{
#' \item H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York, 2009.
#' \item Santander Meteorology Group (2015). downscaleR: Climate data manipulation and
#' statistical downscaling. R package version 0.6-0.
#' https://github.com/SantanderMetGroup/downscaleR/wiki
#' }
#' 
#' 
#' 

getEnsem_comb <- function(..., list = NULL, nrow = 1, legend = TRUE, x = '', y = '', title = '', 
                          output = FALSE) {
  
  if (!is.null(list)) {
    checkBind(list, 'rbind')
    data_ggplot <- do.call('rbind', list)
  } else {
    plots <- list(...)
    checkBind(plots, 'rbind')
    data_ggplot <- do.call('rbind', plots)
  }  
  #data_ggplot$name <- factor(data_ggplot$name, levels = data_ggplot$name, ordered = TRUE)
  
  if (!class(data_ggplot) == 'data.frame') {
    warning('Your input is probably a list, but you forget to add "list = " before it.
            Try again, or check help for more information.')
  } else if (is.null(data_ggplot$name)) {
    stop('No "Name" column in the input data, check the arguments in getFreEnsem() or getHisEnsem(), if 
         output = "ggplot" is assigned, more info please check ?getFreEnsem() or ?getHisEnsem().')
  }
  
  colors = c('brown1', 'dodgerblue3', rainbow(n = length(unique(data_ggplot$variable)) - 2,
                                              start = 0.1))
  
  theme_set(theme_bw())
  mainLayer <- with(data_ggplot, {
    ggplot(data = data_ggplot) +
      aes(x = Date, y = value, color = variable) +
      geom_line(size = 0.5) +
      geom_line(data = data_ggplot[data_ggplot$variable == 'Mean', ], size = 1.6) +
      geom_line(data = data_ggplot[data_ggplot$variable == 'Observation', ], size = 1.6) +
      geom_point(data = data_ggplot[data_ggplot$nav == 1, ], size = 2, shape = 4, color = 'black') +
      scale_color_manual(values = colors) +
      theme(axis.text.x = element_text(size = rel(1.8)),
            axis.text.y = element_text(size = rel(1.8)),
            axis.title.x = element_text(size = rel(1.8)),
            axis.title.y = element_text(size = rel(1.8))) +
      facet_wrap( ~ name, nrow = nrow) +
      labs(x = x, y = y, title = title)
    
  })
  if (legend == FALSE) {
    mainLayer <- mainLayer + 
      theme(legend.position = 'none')
# following ones are to add label, may be added in future.
#      geom_text(data = data_ggplot[data_ggplot$Date == '2003-12-10', ], aes(label = variable), hjust = 0.7, vjust = 1)
#      geom_text(data = data_ggplot[data_ggplot$variable == 'Mean', ], aes(label = variable), hjust = 0.7, vjust = 1)
  }
  
  
  print(mainLayer)
  
  if (output == TRUE) return(data_ggplot)
  
}