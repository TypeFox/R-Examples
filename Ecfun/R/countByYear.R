countByYear <- function(start, end, total=1){
##
## 1.  Confirm class Date 
##
  if(!inherits(start, 'Date')){
    stop('start must be a Date;  has class = ', 
         class(start))
  }
  if(!inherits(end, 'Date')){
    stop('end must be a Date;  has class = ', 
         class(end))
  }
  if(start>end){
    stop('end = ', end, ' before start = ', start)
  }
##
## 2.  compute start, end year 
##
  startYr <- as.integer(substring(start, 1, 4))
  endYr <- as.integer(substring(end, 1, 4))
##
## 2.  startYr==endYr
##
  if(startYr==endYr){
    out <- total
    names(out) <- startYr
  } else {
##
## 3.  days by year 
##
    yrs1 <- startYr:(endYr+1)
    Jan1 <- as.Date(paste(yrs1, '01-01', sep='-'))
    yrDays <- diff(Jan1)
    yrDays[1] <- (Jan1[2]-start)
    nyrs <- length(yrDays)
    yrDays[nyrs] <- (end-Jan1[nyrs]+1)
##
## 4.  allocate
##
    totDays <- as.numeric(sum(yrDays))
    out <- (total * as.numeric(yrDays)/totDays)
    names(out) <- startYr:endYr 
  } 
##
## 5.  Done 
##
  out
}
