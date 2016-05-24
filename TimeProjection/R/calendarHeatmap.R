monthyear = function(timeprojection) {
    paste(timeprojection$month, timeprojection$year)
}

#' Calendar Heatmap
#'
#' Create a plot mimicing a calendar with a heatmap of values
#' @param dates a vector of date objects
#' @param values a numeric vector with same length as dates
#' @examples
#'    dates = timeSequence(from = '2012-01-01', to = '2012-12-31', by = 'day')
#'    plotCalendarHeatmap(as.Date(dates), 1:366)
#' @export
plotCalendarHeatmap = function(dates, values) {
    if(!require (ggplot2) & !require (plyr))
        stop("The packages ggplot2 and plyr are required to use plotCalendarHeatmap")
    tp = projectDate(dates, drop = F)
    tp$values = values
    tp$week = as.numeric(format(dates, "%W"))
    tp = ddply(tp, .(year, month), transform, monthweek=1+week-min(week))
  
    ggplot(tp, aes(monthweek, weekday, fill = values)) + 
        geom_tile(colour = "white") + facet_grid(year~month) +
        scale_fill_gradientn(colours = c("#D61818","#FFAE63","#FFFFBD","#B5E384"))
}

