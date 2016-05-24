#' Get mean rainfall data.
#' 
#' Get mean rainfall data, e.g. mean annual rainfall, mean monthly rainfall and mean winter rainfall.
#' 
#' @param inputTS A time series with only data column (1 column).
#' @param method A string showing the method used to calculate mean value, e.g., "annual".
#' more information please refer to details.
#' @param yearIndex A NUMERIC ARRAY showing the year index of the time series.
#' @param monthIndex A NUMERIC ARRAY showing the month index of the time series.
#' @param fullResults A boolean showing whether the full results are shown, default is FALSE. If 
#' FALSE, only mean value will be returned, if TRUE, the sequence of values will be returned.
#' @param omitNA A boolean showing in the calculation, whether NA is omitted, default is FALSE.
#' @param plot A boolean showing whether the results will be plotted.
#' @param ..., \code{title, x, y} showing the title and x and y axis of the plot, shoud be a string.
#' @details
#' There are following methods to be selected, 
#' "annual": annual rainfall of each year is plotted.  
#' "winter", "spring", "autumn", "summer": seasonal rainfall of each year is plotted.
#' Month(number 1 to 12): month rainfall of each year is plotted, e.g. march rainfall of each year.
#' "meanMonthly": the mean monthly rainfall of each month over the whole period.
#' 
#' Since "winter" is a crossing year, 12, 1, 2, 12 is in former year, and 1, 2 are in latter year.
#' so winter belongs to the latter year.
#' 
#' @return The mean value of the input time series or the full results before calculating mean.
#' 
# data(testdl)
# TS  <- testdl[[1]]
# year = as.numeric(format(TS[, 1], '%Y'))
# month = as.numeric(format(TS[, 1], '%m'))
# 
# # Get the mean spring precipitation.
# a <- getMeanPreci(TS[, 2], method = 'spring', yearIndex = year, monthIndex = month)
# a
# 
# # Get the series of spring precipitation, set fullResults = TRUE.
# a <- getMeanPreci(TS[, 2], method = 'spring', yearIndex = year, monthIndex = month,
#                   fullResults = TRUE)
# a
# 
# # If missing value is excluded, set omitNA = TRUE.
# a <- getMeanPreci(TS[, 2], method = 'winter', yearIndex = year, monthIndex = month,
#                   omitNA = TRUE, fullResults = TRUE)
# a
# 
# # Get special month precipitation, e.g. march.
# a <- getMeanPreci(TS[, 2], method = 3, yearIndex = year, monthIndex = month,
#                   fullResults = TRUE)
# a
# 
# # We can also get annual precipitation.
# a <- getMeanPreci(TS[, 2], method = 'annual', yearIndex = year, monthIndex = month,
#                   fullResults = TRUE)

getMeanPreci <- function(inputTS, method = NULL, yearIndex = NULL, monthIndex = NULL,
                         fullResults = FALSE, omitNA = TRUE, plot = FALSE, ...) {
  # First check if all the records are NA.
  if (any(!is.na(inputTS))) {
    #converting daily preci to the wanted preci.
    if (method == 'annual') {
      ###yearIndex <- startTime$year + 1900
      annualPreci <- tapply(inputTS, INDEX = yearIndex, FUN = sum, na.rm = omitNA)#ggplot is able not to show NA, so choose TRUE
      if (fullResults == TRUE) output <- annualPreci else output <- mean(annualPreci, na.rm = TRUE)
      
    } else if (method == 'meanMonthly') {
      
      monthlypreci <- tapply(inputTS, INDEX = list(yearIndex, monthIndex), FUN = sum, na.rm = omitNA)
      meanMonthlyPreci <- apply(monthlypreci, MARGIN = 2, FUN = mean, na.rm = TRUE)
      
      if (fullResults == TRUE) output <- meanMonthlyPreci else output <- mean(meanMonthlyPreci, na.rm = TRUE)
      
    }else if (method == 'winter') {
#       #winter is the most tricky part, because it starts from Dec to Feb next year, it's a year-crossing season,
#       #so we have to make some changes to the monthIndex
#       #e.g.data from 1950.1.1 - 2008.3.31 if we want to calculate the mean winter preci, to calculate winter month
#       #December, we have to move the yearIndex one month forwards or two months backwards, to make 12,1,2 in one year      
#       ###yearIndex <- startTime$year + 1900
#       ###monthIndex <- startTime$mon + 1
#       
#       #we move the yearIndex one month backwards
#       yearIndex_new <- c(yearIndex[32:length(yearIndex)], rep(tail(yearIndex, 1), 31))
#       
#       winterIndex <- which(monthIndex == 12 | monthIndex == 1 | monthIndex == 2)
#       winterYear <- yearIndex_new[winterIndex]#this index is used for calculation
#       
#       #because we don't have 1949.Dec, so the first winter is not intact, so first two months are elemenated
#       
#       startIndex <- length(which(winterYear == yearIndex[1])) + 1
#       winterOfLastYear <- length(which(winterYear == tail(yearIndex, 1)))
#       if (winterOfLastYear > 91) {
#         endIndex <- length(winterYear) - 31 #in case the data set ends at Dec.31
#         
#       } else if (winterOfLastYear < 90) { # incase the data ends at Jan 31
#         endIndex <- length(winterYear) - length(which(winterYear == tail(yearIndex, 1)))
#         
#       } else {
#         endIndex <- length(winterYear)
#       }
#       
#       inputTS <- inputTS[winterIndex][startIndex:endIndex]#needs two process with inputPreci, first, extract
#       #the winter preci, second, delete first two month of 1950
#       
#       winterYear <- winterYear[startIndex:endIndex]#needs one process, delete two months   
#       seasonalPreci <- tapply(inputTS,INDEX = winterYear, FUN = sum, na.rm = omitNA)
      
      # upper part is the older method saved as backup.
      
      matrix <- tapply(inputTS, INDEX = list(monthIndex, yearIndex), FUN = sum, na.rm = omitNA)
      col <- colnames(matrix)
      dec <- matrix['12',] # extract December.
      dec <- c(NA, dec[1:length(dec) - 1]) # rearrange December order to push it to next year.
      names(dec) <- col
      matrix <- rbind(dec, matrix[rownames(matrix) != '12', ])      
      seasonalPreci <- apply(matrix, MARGIN = 2, function(x) sum(x[c('dec', '1', '2')]))
      
      if (fullResults == TRUE) output <- seasonalPreci else output <- mean(seasonalPreci, na.rm = TRUE)  
      
    } else if (method == 'spring') {
      
#       springIndex <- which(monthIndex == 3 | monthIndex == 4 | monthIndex == 5)
#       springYear <- yearIndex[springIndex]
#       inputTS <- inputTS[springIndex]
#       seasonalPreci <- tapply(inputTS, INDEX = springYear, FUN = sum, na.rm = omitNA)
      
      
      matrix <- tapply(inputTS, INDEX = list(monthIndex, yearIndex), FUN = sum, na.rm = omitNA)
      seasonalPreci <- apply(matrix, MARGIN = 2, function(x) sum(x[c('3', '4', '5')]))
      
      if (fullResults == TRUE) output <- seasonalPreci else output <- mean(seasonalPreci, na.rm = TRUE)
      
    } else if (method == 'summer') {
      
      matrix <- tapply(inputTS, INDEX = list(monthIndex, yearIndex), FUN = sum, na.rm = omitNA)
      seasonalPreci <- apply(matrix, MARGIN = 2, function(x) sum(x[c('6', '7', '8')]))
      
      if (fullResults == TRUE) output <- seasonalPreci else output <- mean(seasonalPreci, na.rm = TRUE)
      
    } else if (method == 'autumn') {
      
      matrix <- tapply(inputTS, INDEX = list(monthIndex, yearIndex), FUN = sum, na.rm = omitNA)
      seasonalPreci <- apply(matrix, MARGIN = 2, function(x) sum(x[c('9', '10', '11')]))

      if (fullResults == TRUE) output <- seasonalPreci else output <- mean(seasonalPreci, na.rm = TRUE)
      
    } else if (is.numeric(method)) {
      
      month <- method
      
      #check if month exist 
      e <- match(month, unique(monthIndex))
      if (is.na(e)) {
        e1 <- paste(unique(monthIndex), collapse = ',')
        m <- paste('No input month exists in the dataset, choose among', e1)
        stop(m)
      }
      
      monthlyPreci <- tapply(inputTS, INDEX = list(yearIndex, monthIndex), 
                             FUN = sum, na.rm = omitNA)[, toString(month)]
      
      if (fullResults == TRUE) output <- monthlyPreci else output <- mean(monthlyPreci, na.rm = TRUE)
    }
    
  } else {
    output <- NA
  }

  if (plot == TRUE) {
    a <- data.frame(Date = names(output), value = output)
    
    theme_set(theme_bw())
    mainLayer <- with(a, {
      ggplot(a) +
        geom_bar(aes(x = Date, y = value), stat = 'identity', fill = 'cyan') +
        labs(empty = NULL, ...) +#in order to pass "...", arguments shouldn't be empty.
        theme(plot.title = element_text(size = rel(1.3), face = 'bold'),
              axis.title.x = element_text(size = rel(1.2)),
              axis.title.y = element_text(size = rel(1.2))) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    })
      
    print (mainLayer)
  }
  
  return(output)
}
