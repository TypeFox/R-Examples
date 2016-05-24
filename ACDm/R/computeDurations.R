computeDurations <- function(transactions, open = "10:00:00", close = "18:25:00", rm0dur = TRUE, type = "trade", priceDiff = .1, cumVol = 10000){ 
  
  open <- as.POSIXlt(strptime(open, "%H:%M:%S"))
  open <- open$h * 3600 + open$min * 60 + open$sec  
  close <- as.POSIXlt(strptime(close, "%H:%M:%S"))
  close <- close$h * 3600 + close$min * 60 + close$sec
  
  type <- switch(type, trade = 1, transactions = 1, price = 2, volume = 3)
  
  
  if("data.frame" %in% class(transactions)){
    
    if(length(transactions$time) == 0) stop("the data.frame 'transactions' must contain a column named 'time' 
                                            with timestamps of each transaction")
    
    if(!("POSIXlt" %in% class(transactions$time))) transactions$time <- as.POSIXlt(transactions$time)
    
  } else{
    
    transactions <- data.frame(time = transactions)
    if(!("POSIXlt" %in% class(transactions$time))) transactions$time <- as.POSIXlt(transactions$time)
    
  }
  

  if(length(transactions$volume) != 0 || length(transactions$price)){ #volume and/or price were provided
    
    temp <- .C("computeDurationsSubSec", 
               as.integer(transactions$time$year), #1
               as.integer(transactions$time$mon),
               as.integer(transactions$time$mday),
               as.integer(transactions$time$hour),
               as.integer(transactions$time$min), #5
               as.double(transactions$time$sec),
               as.integer(rep(0,length(transactions$time))),
               as.integer(rep(0,length(transactions$time))),
               as.integer(rep(0,length(transactions$time))),
               as.integer(rep(0,length(transactions$time))), #10
               as.integer(rep(0,length(transactions$time))),
               as.double(rep(0,length(transactions$time))),
               as.integer(transactions$volume),
               as.double(transactions$price), 
               as.integer(rep(0,length(transactions$time))), #15
               as.double(rep(0,length(transactions$time))),
               as.integer(rep(0,length(transactions$time))),
               as.double(rep(0,length(transactions$time))),
               as.integer(length(transactions$time)),
               as.integer(0), #20
               as.double(open),
               as.double(close),
               as.integer(type),
               as.integer(rm0dur),
               as.double(priceDiff),
               as.integer(cumVol), PACKAGE = "ACDm") #26
    
    
    n <- temp[[20]]
    times <- paste(temp[[7]][1:n] + 1900, temp[[8]][1:n] + 1, temp[[9]][1:n], temp[[10]][1:n], temp[[11]][1:n], temp[[12]][1:n], sep = ":")
    
    dftemp <- data.frame(time =  strptime(times, "%Y:%m:%d:%H:%M:%OS"))
    if(length(transactions$price) != 0) dftemp <- cbind(dftemp,  price = temp[[16]][1:n])
    if(length(transactions$volume) != 0) dftemp <- cbind(dftemp,  volume = temp[[15]][1:n])
    if(rm0dur) dftemp <- cbind(dftemp,  Ntrans = temp[[17]][1:n])
    dftemp <- cbind(dftemp, durations = temp[[18]][1:n])
    
  } else{ #only transaction times were given
    
    temp <- .C("computeDurationsShort", 
               as.integer(transactions$time$year), #1
               as.integer(transactions$time$mon),
               as.integer(transactions$time$mday),
               as.integer(transactions$time$hour),
               as.integer(transactions$time$min), #5
               as.double(transactions$time$sec),
               as.integer(rep(0,length(transactions$time))),
               as.integer(rep(0,length(transactions$time))),
               as.integer(rep(0,length(transactions$time))),
               as.integer(rep(0,length(transactions$time))), #10
               as.integer(rep(0,length(transactions$time))),
               as.double(rep(0,length(transactions$time))),
               as.double(rep(0,length(transactions$time))),     
               as.integer(0), 
               as.integer(rep(0,length(transactions$time))),    #15   
               as.integer(length(transactions$time)),
               as.integer(open),
               as.integer(close),
               as.integer(rm0dur), PACKAGE = "ACDm")  #19
    
    n <- temp[[14]] 
    
    times <- paste(temp[[7]][1:n]+1900, temp[[8]][1:n]+1, temp[[9]][1:n], temp[[10]][1:n], temp[[11]][1:n], temp[[12]][1:n], sep = ":")
    
    dftemp <- data.frame(time =  strptime(times, "%Y:%m:%d:%H:%M:%OS"))
    if(rm0dur) dftemp <- cbind(dftemp,  Ntrans = temp[[15]][1:n])
    dftemp <- cbind(dftemp, durations = temp[[13]][1:n])
    
  }
  
  #checks if any of the durations are negative:
  if(any(transactions$durations < 0)){
    
    if(is.unsorted(transactions$time)){
      warning("the provided 'time' column is not in chronological order")
    } else{
      warning("Negative durations computed.")
    }
    
  } 
  
  cat("The", length(transactions$time), "transactions resulted in", n, "durations")
  return(dftemp)
  
}