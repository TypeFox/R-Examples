countsByYear <- function(data, start="Start1", 
      end='End1', total='BatDeath', event='WarName', 
      endNA=max(data[, c(start,end)])){
##
## 1.  checks 
##
  if(!inherits(data, 'data.frame')){
    stop('data must be a data.frame;  has class = ', 
         class(data))
  }
  ns <- sum(start == names(data))
  if(ns<1){
    stop('start = ', start, ' not found in names(data)')
  }
  if(ns>1){
    stop('start = ', start, ' found ', ns,
         ' times in names(data)')
  }
  ne <- sum(end == names(data))
  if(ne<1){
    stop('end = ', end, ' not found in names(data)')
  }
  if(ne>1){
    stop('end = ', end, ' found ', ne, ' times in names(data)')
  }
  nv <- (sum(event == names(data)))
  if(nv<1){
    stop('event = ', event, ' not found in names(data)')
  }
  if(nv>1){
    stop('event = ', event, ' found ', ne, ' times in names(data)')
  }
  nr <- nrow(data)
#
  if(!inherits(data[, start], 'Date')){
    stop('data[, start] must be a Date;  has class = ', 
       class(data[, start]))
  }
  if(!inherits(data[, end], 'Date')){
    stop('data[, end] must be a Date;  has class = ', 
         class(data[, end]))
  }
  if(any(is.na(data[, start]))){
      stop('NA found in start;  not allowed')    
  }
  endNAs <- is.na(data[, end])
  data[endNAs, end] <- endNA
##
## 2.  compute start, end year 
##
  firstStart <- min(data[, start])
  startYr <- as.integer(substring(firstStart, 1, 4))
  lastEnd <- max(data[, end])
  endYr <- as.integer(substring(lastEnd, 1, 4))
  yrs <- startYr:endYr 
  ny <- length(yrs) 
  Counts <- matrix(0, ny, nr, 
        dimnames=list(yrs, data[, event]))
##
## 3.  loop by event 
##
  for(i in 1:nr){
    counti <- countByYear(data[i, start], data[i, end], 
                          data[i, total])  
    Counts[names(counti), i] <- as.numeric(counti)
  }
##
## 5.  Done 
##
  Counts
}
