#' This function identifies chart signals in an [0,100]-Points Bullish Percent Chart
#'
#' @param data Input data
bp.signalprocessor <- function(data) {
  for (i in 1:nrow(data)) {
    if (data$status.xo[i]=="X") {
      # we are in x-mode
      if (data$status.bs[i]=="Buy") {
        # we have a buy signal
        data$signal[i] <- "Bull Confirmed"
      } else {
        # we have a sell signal
        if (data$high[i] <= 30 | (i>1 & data$signal[i-1]=="Bull Alert")) {
          data$signal[i] <- "Bull Alert"
        } else {
          data$signal[i] <- "Bear Correction"
        }
      }
    } else {
      # we are in O-mode
      if (data$status.bs[i]=="Sell") {
        # we have a sell signal
        data$signal[i] <- "Bear Confirmed"
      } else {
        # we have a buy signal
        if (data$high[i] >= 70 | (i>1 & data$signal[i-1]=="Bear Alert")) {
          data$signal[i] <- "Bear Alert"
        } else {
          data$signal[i] <- "Bull Correction"
        }
      }
    }
  }
  data
}
