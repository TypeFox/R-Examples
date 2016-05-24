#' plot time series, with marks on missing value.
#' 
#' @param ... input time series.
#' @param type A string representing the type of the time series, e.g. 'line' or 'bar'.
#' @param output A string showing which type of output you want. Default is "data", if "ggplot", the 
#' data that can be directly plotted by ggplot2 will be returned, which is easier for you to make series
#' plots afterwards. 
#' @param name If \code{output = 'ggplot'}, name has to be assigned to your output, in order to differentiate
#' different outputs in the later multiplot using \code{plotTS_comb}.
#' @param plot representing the plot type, there are two types, "norm" and "cum", "norm" gives an normal
#' plot, and "cum" gives a cumulative plot. Default is "norm".
#' @param x label for x axis.
#' @param y label for y axis.
#' @param title plot title.
#' @param list If your input is a list of time series, then use \code{list = your time sereis list}
#' @return A plot of the input time series.
#' @details 
#' If your input has more than one time series, the program will only plot the common period of 
#' different time series.
#' @examples
#' plotTS(testdl[[1]])
#' plotTS(testdl[[1]], x = 'xxx', y = 'yyy', title = 'aaa')
#' 
#' # If input is a datalist
#' plotTS(list = testdl)
#' 
#' # Or if you want to input time series one by one
#' # If plot = 'cum' then cumulative curve will  be plotted.
#' plotTS(testdl[[1]], testdl[[2]], plot = 'cum')
#' 
#' # You can also directly plot multicolumn dataframe
#' dataframe <- list2Dataframe(extractPeriod(testdl, commonPeriod = TRUE))
#' plotTS(dataframe, plot = 'cum')
#' 
#' # Sometimes you may want to process the dataframe and compare with the original one
#' dataframe1 <- dataframe
#' dataframe1[, 2:4] <- dataframe1[, 2:4] + 3
#' plotTS(dataframe, dataframe1, plot = 'cum')
#' # But note, if your input is a multi column dataframe, it's better to plot one using plotTS,
#' # and compare them using plotTS_comb. If all data are in one plot, there might be too messy.
#' 
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @references 
#' \itemize{
#' \item H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York, 2009.
#' }
#' 
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
plotTS <- function(..., type = 'line', output = 'data', plot = 'norm', name = NULL, x = NULL, 
                   y = NULL, title = NULL, list = NULL) {
  ## arrange input TS or TS list.
  if (is.null(list)) {
    list <- list(...)
    if (!class(list[[1]]) == 'data.frame') {
      warning('Your input is probably a list, but you forget to add "list = " before it.
              Try again, or check help for more information.')
    }
#     Following part is for plot different time series with different date, but too complicated
#     using ggplot. and normal use doesn't need such process. So save it as backup.
#     listNames <- names(list)
#     # in order to be used later to differentiate lists, there should be a name for each element.
#     # Then assign the name column to each list element.
#     if (is.null(listNames)) listNames <- 1:length(list)
#     
#     giveName <- function(x, y) {
#       colnames(x) <- NULL
#       x$TSname <- rep(listNames[y], nrow(x))
#       return(x)
#     }
#     list1 <- mapply(FUN = giveName, x = list, y = 1:length(list), SIMPLIFY = FALSE)
#     
#     checkBind(list1, 'rbind')
#     
#     TS <- do.call('rbind', list1)
  }
  
  list_common <- extractPeriod(list, commonPeriod = TRUE)
  TS <- list2Dataframe(list_common)
  
  if (!is.null(names(list)) & (ncol(TS) - 1) == length(list)) colnames(TS)[2:(length(list) + 1)] <- names(list)
  
  # Check input, only check the first column and first row.
  if (!grepl('-|/', TS[1, 1])) {
    stop('First column is not date or Wrong Date formate, check the format in ?as.Date{base} 
         and use as.Date to convert.')
  }
  
  TS[, 1] <- as.Date(TS[, 1])
  colnames(TS)[1] <- 'Date'
  
  # first column's name may not be Date, so change its name to Date
  
  data_plot <- melt(TS, id.var = 'Date')
  NAIndex <- which(is.na(data_plot$value))
  
  # assign 0 to NA values
  if (plot == 'norm') {
    data_plot$value[NAIndex] <- 0
    lineSize <- 0.7
  } else if (plot == 'cum') {
    TS[is.na(TS)] <- 0
    cum <- cbind(data.frame(Date = TS[, 1]), cumsum(TS[2:ncol(TS)]))
    
    data_plot <- melt(cum, id.var = 'Date')
    lineSize <- 1
  }
  
  
  # Assigning x, y and title
  if (is.null(x)) x <- colnames(TS)[1]
  # y aixs cannot decide if it's a multi column dataframe
  #if (is.null(y)) y <- names[2]
  
  theme_set(theme_bw())
  mainLayer <- with(data_plot, {
    ggplot(data = data_plot) +
    # It's always better to use colname to refer to
      aes(x = Date, y = value, color = variable) +
      theme(plot.title = element_text(size = rel(1.8), face = 'bold'),
            axis.text.x = element_text(size = rel(1.8)),
            axis.text.y = element_text(size = rel(1.8)),
            axis.title.x = element_text(size = rel(1.8)),
            axis.title.y = element_text(size = rel(1.8))) +
      labs(x = x, y = y, title = title)
  })
  
  
#  color <- 'dodgerblue4'
  if (type == 'bar') {
    secondLayer <- with(data_plot, {
      geom_bar(stat = 'identity')
    })
  } else if (type == 'line') {
    secondLayer <- with(data_plot, {
      geom_line(size = lineSize)
    })
  } else {
    stop("No such plot type.")
  }
  
  
  missingVLayer <- with(TS, {
    geom_point(data = data_plot[NAIndex, ], group = 1, size = 3, shape = 4, color = 'black')
  })
  
  plotLayer <- mainLayer + secondLayer + missingVLayer
  
  print(plotLayer) 
  
  if (output == 'ggplot') {
    if (is.null(name)) stop('"name" argument not found, 
                            If you choose "ggplot" as output, please assign a name.')
    
    data_plot$name <- rep(name, nrow(data_plot))     
    data_plot$nav <- rep(0, nrow(data_plot))
    data_plot$nav[NAIndex] <- 1
    return(data_plot)
  }
}




#' Combine time seires plot together
#' @param ... different time series plots generated by \code{plotTS(, output = 'ggplot')}, refer to details.
#' @details
#' ..., representing different ouput file generated by \code{plotTS(, output = 'ggplot'), name = yourname}, 
#' different names must be assigned when generating different output.
#' 
#' e.g.
#' a1, a2, a3 are different files generated by \code{plotTS(, output = 'ggplot'), name = yourname}, you can
#' set \code{plotTS(a1,a2,a3)} or \code{plotTS(list = list(a1,a2,a3))}
#' 
#' @param nrow A number showing the number of rows.
#' @param type A string showing 'line' or 'bar'.
#' @param list If input is a list containing different ggplot data, use l\code{list = inputlist}.
#' @param x A string of x axis name.
#' @param y A string of y axis name.
#' @param title A string of the title.
#' @param output A boolean, if chosen TRUE, the output will be given.
#' NOTE: yOU HAVE TO PUT A \code{list = }, before your list.
#' @return A combined time series plot.
#' @examples
#' a1 <- plotTS(testdl[[1]], output = 'ggplot', name = 1)
#' a2 <- plotTS(testdl[[2]], output = 'ggplot', name = 2)
#' 
#' plotTS_comb(a1, a2)
#' plotTS_comb(list = list(a1, a2), y = 'y axis', nrow = 2)
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @references 
#' \itemize{
#' \item H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York, 2009.
#' }
#' @export
#' @import ggplot2
plotTS_comb <- function(..., nrow = 1, type = 'line', list = NULL, x = 'Date', y = '', title = '', 
                        output = FALSE){
  # In ploting the time series, since the data comes from outside of hyfo, 
  # It's more complicated, since they don't always have the same
  # column name, if not, there is not possible to do rbind.
  # So we need to first save the name, and rbind, and put back the name.
  
  if (!is.null(list)) {
    checkBind(list, 'rbind')
    data_ggplot <- do.call('rbind', list)
  } else {
    
    bars <- list(...)
    checkBind(bars, 'rbind')
    data_ggplot <- do.call('rbind', bars)
  }
  
  if (!class(data_ggplot) == 'data.frame') {
    warning('Your input is probably a list, but you forget to add "list = " before it.
            Try again, or check help for more information.')
  } else if (is.null(data_ggplot$name)) {
    stop('No "name" column in the input data, check the arguments in getPreciBar(), if 
         output = "ggplot" is assigned, more info please check ?getPreciBar.')
  }

  
  theme_set(theme_bw())
  mainLayer <- with(data_ggplot, {
    ggplot(data = data_ggplot) +
      # It's always better to use colname to refer to
      aes(x = Date, y = value, color = variable) +
      theme(plot.title = element_text(size = rel(1.8), face = 'bold'),
            axis.text.x = element_text(angle = 90, hjust = 1, size = rel(1.8)),
            axis.text.y = element_text(size = rel(1.8)),
            axis.title.x = element_text(size = rel(1.8)),
            axis.title.y = element_text(size = rel(1.8))) +
      geom_point(data = data_ggplot[data_ggplot$nav == 1, ], size = 2, shape = 4, color = 'red') +
      facet_wrap( ~ name, nrow = nrow) +
      labs(x = x, y = y, title = title)
    
  })
  
  
  if (type == 'bar') {
    secondLayer <- with(data_ggplot, {
      geom_bar(stat = 'identity', size = 1)
    })
  } else if (type == 'line') {
    secondLayer <- with(data_ggplot, {
      geom_line(size = 1)
    })
  } else {
    stop("No such plot type.")
  }
  
  print(mainLayer + secondLayer)
  
  if (output == TRUE) return(data_ggplot)
}




#' get L moment analysis of the input distribution
#' 
#' @param dis A distribution, for hydrology usually a time series with only data column without time.
#' @return The mean, L-variation, L-skewness and L-kurtosis of the input distribution
#' @examples
#' dis <- seq(1, 100)
#' getLMom(dis)
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @export
#' @references 
#' 
#' \itemize{
#' \item J. R. M. Hosking (2015). L-moments. R package, version 2.5. URL:
#' http://CRAN.R-project.org/package=lmom.
#' }
#' 
#' 
#' @importFrom lmom samlmu
#' 
getLMom <- function(dis){
  
  LMom <- samlmu(dis, nmom = 4, ratios = TRUE)
  
  mean <- LMom[1]
  LCV <- LMom[2]/LMom[1]
  Lskew <- LMom[3]
  Lkur <- LMom[4]
  
  output <- data.frame(mean = mean, Lcv = LCV, Lskew = Lskew, Lkur = Lkur)
  return(output)
}

#' get moment analysis of the input distribution
#' 
#' @param dis A distribution, for hydrology usually a time series with only data column without time.
#' @return The mean, variation, skewness and kurtosis of the input distribution
#' @examples
#' dis <- seq(1, 100)
#' getMoment(dis)
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @export
#' @references 
#' 
#' \itemize{
#' \item Lukasz Komsta and Frederick Novomestky (2015). moments: Moments, cumulants, skewness, kurtosis and
#' related tests. R package version 0.14. http://CRAN.R-project.org/package=moments
#' 
#' \item R Core Team (2015). R: A language and environment for statistical computing. R Foundation for
#' Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.
#' }
#' 
#' @importFrom moments skewness kurtosis
#' @importFrom stats var
getMoment <- function(dis) {
  mean <- mean(dis, na.rm = TRUE)
  variance <- var(dis, na.rm = TRUE)
  skewness <- skewness(dis, na.rm = TRUE)
  kurtosis <- kurtosis(dis, na.rm = TRUE)
  
  output <- data.frame(mean=mean, Variance = variance, Skewness = skewness, Kurtosis = kurtosis)
  
  return(output)
}
