# This function analyzes a (preliminary) P&F Chart for Bullish Support Line and Bearish Resistance Line
# 
# Finding the appropriate trendlines is explained very good at \url{http://stockcharts.com/school/doku.php?id=chart_school:chart_analysis:pnf_charts:pnf_trendlines}.
# 
# @seealso \url{http://stockcharts.com/school/doku.php?id=chart_school:chart_analysis:pnf_charts:pnf_trendlines}
xo.trendline.processor <- function(data, slope=1) {
  # define local function which determines the slope of the trend lines
  trend.offset <- function(start.column,recent.column) {
    floor((recent.column-start.column)*slope)
  }
  data$tl.brl.boxnumber <- rep(x=NA,times=nrow(data)) # initialize bearish resistance line 
  data$tl.bsl.boxnumber <- rep(x=NA,times=nrow(data)) # initialize bullsih support line 
  data$tl.status <- rep(x=NA,times=nrow(data)) # initialize trendline indicator
  # as all charts start with X-column by default, we start with a bearish resistance line
  trend.start.column <- min(data$column,na.rm=T) # column number, where trendline starts
  trend.start.boxnumber <- max(data$boxnumber[data$column <= trend.start.column],na.rm=T)+1 # boxnumber, where current trend starts
  trend.status.bullish <- FALSE # current trendline is bullish
  
  for (i in 1:nrow(data)) {
    # processing may start erliest at second column
    if (data$column[i]>min(data$column)) {
      if (trend.status.bullish) {
        trend.current.boxnumber <- trend.start.boxnumber + trend.offset(trend.start.column,data$column[i])
        # check if current trendline is still valid, else start new bearish trend
        if (data$boxnumber[i] >= trend.current.boxnumber) {
          data$tl.bsl.boxnumber[i] <- trend.current.boxnumber 
          data$tl.status[i] <- "BULLISH"
        } else {
          # trend was broken, create new bearish trend
          trend.status.bullish <- FALSE
          trend.start.column <- data$column[i]
          trend.start.boxnumber <- max(data$boxnumber[data$column == trend.start.column-1],na.rm=T)
          data$tl.brl.boxnumber[i] <- trend.start.boxnumber
          data$tl.status[i] <- "BEARISH"
        }
      } else {
        trend.current.boxnumber <- trend.start.boxnumber - trend.offset(trend.start.column,data$column[i])
        # check if current trendline is still valid, else start new bullish trend
        if (data$boxnumber[i] <= trend.current.boxnumber) {
          data$tl.brl.boxnumber[i] <- trend.current.boxnumber
          data$tl.status[i] <- "BEARISH"
        } else {
          # trend was broken, create new bullish trend
          trend.status.bullish <- TRUE
          trend.start.column <- data$column[i]
          trend.start.boxnumber <- min(data$boxnumber[data$column == trend.start.column-1],na.rm=T)
          data$tl.bsl.boxnumber[i] <- trend.start.boxnumber
          data$tl.status[i] <- "BULLISH"
        }
      }
    }
  }
  # return result
  data
}
