#' Draw an interactive line plot for cumulative P&L 
#'
#' This function takes a data frame and returns an interactive plot for the
#' cumulative P&L of the portfolio data. The plot is zoomable and specific data can be shown once the 
#' cursor is placed on the graph. 
#' 
#' @param x A data frame that the function takes in. 
#' 
#' @return A plot. Red bars indicate loss and green bars indicate profit.

cumpnl_plot <- function(x) {
  
  ## Fool the R CMD check
  
  date <- name <- cumpnl <- NULL
  
  title.var <- paste("Cumulative P&L for", x$name[1])
  
  ## Convert the data frame to a time-series class so that dygraphs can accept
  ## the data frame and generate a chart
  
  x <- xts(x$cumpnl, order.by = x$date)
  
  ## Do the actual plotting by specifying the data frame and title in "dygraph",
  ## specifying legend formats in "dySeries", add a navigation bar for zooming
  ## in "dyRangeSelector", modify the legend's layout in "dyLegend" and format
  ## the display of numbers in "dyOptions" so that 6,000,000 will be displayed
  ## as 6M.
  
  dygraph(x, main = title.var) %>%
    dySeries(label = "Cumulative P&L") %>%
    dyRangeSelector() %>%
    dyLegend(labelsSeparateLines = TRUE) %>%
    dyOptions(labelsKMB = TRUE,
              useDataTimezone = TRUE,
              fillGraph = TRUE)
}