#' Identifiy for a given P&F Table the current vertical price objective 
#' triggered by the last signal reversal.
#' 
#' @param data Input data
#' @param reversal Number of boxes for reversal
#' @param boxsize Size of one box
#' @param log Use logarithmic scale
currentVPOBreakoutMethod <- function(data,reversal,boxsize,log) {
  price.objective <- list(boxnumber=NA,price=NA)
  if (nrow(data)>=2) {
    # select only reversal days in data
    # FIXME the next statement is the bottleneck in this function!!!
    index <- c(FALSE,data$status.bs[1:(nrow(data)-1)]!=data$status.bs[2:nrow(data)])
    mydata <- data[index,]
    if (nrow(mydata)>0) {
      # identify current status.bs
      current.status <- data$status.bs[nrow(data)]
      # find latest reversal column
      reversal.column <- max(mydata$column,na.rm=T)
      # find minimum and maximum boxnumber of price objective column in original data
      min.boxnumber <- NA
      max.boxnumber <- NA
      if (current.status=="Buy") {
        max.boxnumber <- max(data$boxnumber[data$column==reversal.column],na.rm=T)
        min.boxnumber <- min(data$boxnumber[data$column==reversal.column-1],na.rm=T)+1      
      } else if (current.status=="Sell") {
        min.boxnumber <- min(data$boxnumber[data$column==reversal.column],na.rm=T)
        max.boxnumber <- max(data$boxnumber[data$column==reversal.column-1],na.rm=T)-1        
      } 
      # determine extension estimate
      # extension.estimate.in.boxes <- (max.boxnumber-min.boxnumber)*reversal
      # determine price objective box
      boxnumber <- NA
      price <- NA
      if (current.status=="Buy") {
        boxnumber <- min.boxnumber + (max.boxnumber-min.boxnumber+1)*reversal
        # translate price.objective.box into real number
        price <- box2lower(boxnumber=boxnumber,boxsize=boxsize,log=log)
      } else if (current.status=="Sell") {
        boxnumber <- max.boxnumber - (max.boxnumber-min.boxnumber+1)*(reversal-1)
        # translate price.objective.box into real number
        price <- box2upper(boxnumber=boxnumber,boxsize=boxsize,log=log)
      } else {
        # should not happen
        stop("Internal error in .currentVerticalPriceObjective()!")
      }
      price.objective <- list(boxnumber=boxnumber,price=price)
    }
  } 
  price.objective
}
