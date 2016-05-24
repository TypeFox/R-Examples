#' Analyzes a given PNF time-series for Buy&Sell patterns
#' 
#' @param data Input data
#' @param reversal Number of boxes for reversal
xo.signalprocessor <- function(data, reversal=3) {
  # check for needed columns
  if (!"boxnumber" %in% names(data))
    stop("column 'boxnumber' is missing!")
  if (!"column" %in% names(data))
    stop("column 'column' is missing!")
  if (!"status.xo" %in% names(data))
    stop("column 'status.xo' is missing!")
  if (!"date" %in% names(data))
    stop("column 'date' is missing!")
  if ("signal.bs" %in% names(data))
    warning("column 'signal.bs' already exists, will be overriden!")
  
  #
  # initialize decision tree and first signals
  data$signal.bs[1] <- "UNDEFINED"
  
  #
  # For every row in data go through decision tree
  for (i in 2:nrow(data)){
    if (data$status.xo[i] == "X") {
      # check X-branch of decision tree
      if (raisingTop(data[data$date<=data$date[i],],data$column[i])) {
        # we have at least DOUBLE TOP
        data$signal.bs[i] <- "DOUBLE TOP"
        if (raisingBottom(data[data$date<=data$date[i],],data$column[i]-1)) {
          data$signal.bs[i] <- "BULLISH SIGNAL"
          if (raisingTop(data[data$date<=data$date[i],],data$column[i]-2)) {
            data$signal.bs[i] <- "TRIPLE BULLISH SIGNAL"
            if ((doubleBottom(data[data$date<=data$date[i],],data$column[i]-3) & 
                   doubleTop(data[data$date<=data$date[i],],data$column[i]-4))) {
              data$signal.bs[i] <- "BULLISH CATAPULT"
            } 
          } else if (fallingTop(data[data$date<=data$date[i],],data$column[i]-2)) {
            data$signal.bs[i] <- "BULLISH TRIANGLE"
          }
        } else if (doubleTop(data[data$date<=data$date[i],],data$column[i]-2)) {
          data$signal.bs[i] <- "TRIPLE TOP"
        } else if ((fallingBottom(data[data$date<=data$date[i],],data$column[i]-1) &
                      fallingTop(data[data$date<=data$date[i],],data$column[i]-2) &
                      fallingBottom(data[data$date<=data$date[i],],data$column[i]-3))) {
          data$signal.bs[i] <- "BEARISH SIGNAL REVERSED"
        } 
        # TODO low pole is a tricky one!!
        #       } else if (((data$column[i]>=4) & 
        #                   (minBox(data[data$date<=data$date[i],],data$column[i]-3)-reversal<=minBox(data[data$date<=data$date[i],],data$column[i]-1)) &
        #                   (2*(maxBox(data[data$date<=data$date[i],],data$column[i])-minBox(data[data$date<=data$date[i],],data$column[i]-2))>(maxBox(data[data$date<=data$date[i],],data$column[i]-1)-minBox(data[data$date<=data$date[i],],data$column[i]-1))))) {
        #         data$signal.bs[i] <- "LOW POLE"
      } else if (FALSE) {
        # TODO insert condtions for BEAR TRAP
        data$signal.bs[i] <- "BEAR TRAP"
      } else {
        # default case: previous signal is still valid
        data$signal.bs[i] <- data$signal.bs[i-1]
      }
    } else {
      # check O-branch of decision tree
      if (fallingBottom(data[data$date<=data$date[i],],data$column[i])) {
        # we have at least DOUBLE BOTTOM
        data$signal.bs[i] <- "DOUBLE BOTTOM"
        if (fallingTop(data[data$date<=data$date[i],],data$column[i]-1)) {
          data$signal.bs[i] <- "BEARISH SIGNAL"
          if (fallingBottom(data[data$date<=data$date[i],],data$column[i]-2)) {
            data$signal.bs[i] <- "TRIPLE BEARISH SIGNAL"
            if ((doubleTop(data[data$date<=data$date[i],],data$column[i]-3) &
                   doubleBottom(data[data$date<=data$date[i],],data$column[i]-4))) {
              data$signal.bs[i] <- "BEARISH CATAPULT"
            }
          } else if (raisingBottom(data[data$date<=data$date[i],],data$column[i]-2)) {
            data$signal.bs[i] <- "BEARISH TRIANGLE"
          }
        } else if (doubleBottom(data[data$date<=data$date[i],],data$column[i]-2)) {
          data$signal.bs[i] <- "TRIPLE BOTTOM"
        } else if ((raisingTop(data[data$date<=data$date[i],],data$column[i]-1) & 
                      raisingBottom(data[data$date<=data$date[i],],data$column[i]-2) &
                      raisingTop(data[data$date<=data$date[i],],data$column[i]-3))) {
          data$signal.bs[i] <- "BULLISH SIGNAL REVERSED"
        }
      } else if (FALSE) {
        # TODO insert condtions for HIGH POLE, this is a tricky one!!!
        data$signal.bs[i] <- "HIGH POLE"
      } else if (FALSE) {
        # TODO insert condtions for BULL TRAP
        data$signal.bs[i] <- "BULL TRAP"
      } else {
        # default case: previous signal is still valid
        data$signal.bs[i] <- data$signal.bs[i-1]
      } 
    }
    data$signal.bs[i]
  } # end for (i in 2:nrow(data))
  data
}
