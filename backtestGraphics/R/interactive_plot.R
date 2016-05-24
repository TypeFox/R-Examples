#' Interactive plot with dygraph
#'
#' This function takes a data frame and returns an interactive plot for 
#' the portfolio data. The plot is zoomable and specific data can be shown once the 
#' cursor is placed on the graph. 
#' 
#' For all bar plots, green represents profit and red represents loss. 
#' This function also calculates gmv, cumulative P&L and return rates if 
#' these statistics are not in the data frame. Number of contracts will all be zero if the column 
#' is missing. 
#' 
#' @param x A data frame that the function takes in. 
#' 
#' @param type The value depicted by the plot, choices are \code{nmv}, \code{gmv},
#' \code{contract}, \code{pnl} and \code{ret}.
#' 
#' @return A plot. Red bars indicate loss and green bars indicate profits.

interactive_plot <- function(x, type) {
  
  ## Initiate variables so that CMD Check won't release notes
  
  name <- date <- nmv <- gmv <- pnl <- cumpnl <- contract <- y <- NULL
  
  ## Set up the titles and labels for different graphs based on user's choice 
  ## and store them for graphics. The conditons here first check which data
  ## the function is going to plot. Then, the function assign a full name to
  ## a variable called "full.type" so that later on the function can use this
  ## variable to paste a title. The function also extract the date column and
  ## the target column from the input data frame and assign them with "x" and
  ## "y" as column names, 
  
  ## At the mean time, the "pnl" column is still necessary because the function
  ## has to divide the data frame into two part, depending on whether the pnl
  ## values are positive or negative. For the part where all pnl values are
  ## positive, the corresponding bars in the interactive plot are shown in
  ## green, while the other bars are shown in red, representing negative pnl
  ## on those days.
  
  if ("nmv" %in% type) {
    full.type <- "Net Market Value"
    title.var <- paste(full.type, "for", x$name[1])
    x <- select(x, x = date, y = nmv, pnl = pnl)
  }
  
  if ("gmv" %in% type) {
    full.type <- "Gross Market Value"
    title.var <- paste(full.type, "for", x$name[1])
    x <- select(x, x = date, y = gmv, pnl = pnl)
  }
  
  if ("pnl" %in% type) {
    full.type <- "Profit and Loss"
    title.var <- paste(full.type, "for", x$name[1])
    x <- select(x, x = date, pnl = pnl)
    x <- mutate(x, y = pnl)
  }
  
  if ("cumpnl" %in% type) {
    full.type <- "Cumulative Profit and Loss"
    title.var <- paste(full.type, "for", x$name[1])
    x <- select(x, x = date, y = cumpnl, pnl = pnl)
  }
  
  if ("contract" %in% type) {
    full.type <- "Number of Contracts"
    title.var <- paste(full.type, "for", x$name[1])
    x <- select(x, x = date, y = contract, pnl = pnl)
  }
  
  ## Divide the data frame into two series according to P&L. One series contain
  ## the dates when you are earning profits, and that series is colored in 
  ## green. The other one represents loss and is colored in red.
  
  profit <- mutate(x, y = ifelse(pnl >= 0, y, NA))
  loss   <- mutate(x, y = ifelse(pnl <  0, y, NA))
  
  ## Change the data frames into time-series format so dygraphs can read it. We
  ## perform a column bind to combine the two time series data into one.
  
  graph_data_profit <- xts(x = profit[["y"]], order.by = profit[["x"]])
  
  graph_data_loss   <- xts(x = loss[["y"]],   order.by = loss[["x"]])
  
  graph_data <- cbind(graph_data_loss, graph_data_profit)
  
  ## Tag the two series so that these tags can be shown in the plot legends. When
  ## the user's mouse is hovering above a specific bar, the legend is no longer
  ## "V1: 1.356 M", but instead "Some instrument at profit: 1.356 M"
  
  names(graph_data) <- c(paste(full.type, "at loss"),
                         paste(full.type, "at profit"))
  
  ## The actual plotting function with dygraph, which is fast and interactive.
  ## In the previous developments, we tried different packages. ggplot is good
  ## at plotting and it is handy to assign different colors to different bars,
  ## but this plot is static and it does not look good to use Shiny to
  ## complement the interactive part. Highcharts is good at plotting interactive
  ## plots and we can change the characteristics of each bar by organizing all
  ## the information into a list, but Highcharts is slow in plotting more than a
  ## thousand bars. rCharts package, except for the Highcharts part, is an even
  ## slower package. All other packages except for these three cannot plot many
  ## bars, as far as we explored.
  
  ## 1. "dygraph". We call the dygraph function and throw in the time series
  ## data as well as the title. 
  ## 2. "dyLegend". We define the legends so that the legends are aligned and
  ## does not appear if the user's mouse is outside the plot.
  ## 3. "dyRangeSelector". We ask dygraphs to put down a navigation bar beneath
  ## the plot.
  ## 4. "dyOptions". We define colors for the two columns in the time series
  ## data (remember we did a cbind before?). We also ask dygraph to simplify
  ## the big numbers, so if the original data is "1,356,105", the plot only
  ## shows "1.356 M". The "plotter" thing takes in a clip of Javascript and
  ## directly throw that into the Javascrip kernal of dygraphs to generate bars
  ## in the plot.
  
  p2 <- dygraph(graph_data,
                main = title.var) %>% 
    dyLegend(labelsSeparateLines = TRUE,
             show = "onmouseover") %>%
    dyRangeSelector() %>%
    dyOptions(colors = c("red", "green"),
              labelsKMB = TRUE,
              plotter = 
                "function barChartPlotter(e) {
                var ctx = e.drawingContext;
                var points = e.points;
                var y_bottom = e.dygraph.toDomYCoord(0);  
                
                var bar_width = 1/2 * (points[1].canvasx - points[0].canvasx);
                ctx.fillStyle = e.color;
                
                for (var i = 0; i < points.length; i++) {
                var p = points[i];
                var center_x = p.canvasx;  // center of the bar
                
                ctx.fillRect(center_x - bar_width / 2, p.canvasy,
                bar_width, y_bottom - p.canvasy);
                ctx.strokeRect(center_x - bar_width / 2, p.canvasy,
                bar_width, y_bottom - p.canvasy);
                }
                }")
  
  return(p2)
}