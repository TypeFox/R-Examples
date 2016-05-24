#' Concatenate Time Series and Resolve Overlap Automatically
#' 
#' Append time series to each other. Resolve overlap determines
#'  which of two ts class time series is
#' reaching further and arranges the two series into first and second 
#' series accordingly. Both time series are concatenated to one
#' if both series had the same frequency. Typically this function is used 
#' concatenate two series that have a certain overlap, but one series clearly 
#' starts earlier while the other lasts longer. If one series starts earlier and 
#' stops later, all elements of the shorter series will be inserted into the 
#' larger series, i.e. elements of the smaller series will replace the elements
#' of the longer series. Usually ts2 is kept. 
#'
#' @param ts1 ts time series, typically the older series
#' @param ts2 ts time series, typically the younger series
#' @param keep_ts2 logical should ts2 be kept? Defaults to TRUE.
#' @export
#' @examples
#' ts1 <- ts(rnorm(100),start = c(1990,1),frequency = 4)
#' ts2 <- ts(1:18,start = c(2000,1),frequency = 4)  
#' resolveOverlap(ts1,ts2)   
#' 
#' # automatical detection of correction sequence!
#' ts1 <- ts(rnorm(90),start = c(1990,1),frequency = 4)            
#' ts2 <- ts(1:60,start = c(2000,1),frequency = 4)
#' resolveOverlap(ts1,ts2)   
#' 
#' # both series are of the same length use sequence of arguments.
#' ts1 <- ts(rnorm(100),start = c(1990,1),frequency = 4)
#' ts2 <- ts(1:48,start = c(2003,1),frequency = 4)  
#'resolveOverlap(ts1,ts2)
#' ts1 <- ts(rnorm(101),start = c(1990,1),frequency = 4)            
#' ts2 <- ts(1:61,start = c(2000,1),frequency = 4)
#' resolveOverlap(ts1,ts2)
#' #' clearly dominatn ts2 series
#' ts1 <- ts(rnorm(50),start = c(1990,1),frequency = 4)            
#' ts2 <- ts(1:100,start = c(1990,1),frequency = 4)
#' resolveOverlap(ts1,ts2)
resolveOverlap <- function(ts1,ts2,keep_ts2=T){
  stopifnot(is.ts(ts1))
  stopifnot(is.ts(ts2))
  stopifnot(frequency(ts1) == frequency(ts2))
  freq <- frequency(ts1)
    
  ts1s <- min(time(ts1))
  ts1e <- max(time(ts1))
  
  ts2s <- min(time(ts2))
  ts2e <- max(time(ts2))
  
  # Cases ####
  if(ts1s < ts2s & ts1e < ts2e){
    out <- c(ts1[1:(which(time(ts1) == ts2s)-1)],ts2)
    out <- ts(out,start = ts1s,frequency = freq)
  } else if(ts1s < ts2s & ts1e >= ts2e){
    a <- which(time(ts1) == ts2s)
    b <- which(time(ts1) == ts2e)
    ts1[a:b] <- ts2
    out <- ts1
  } else if(ts1s == ts2s & ts2e >= ts1e){
    out <- ts2
  } else if(ts1s > ts2s & ts1e < ts2e){
    a <- which(time(ts2) == ts1s)
    b <- which(time(ts2) == ts1e)
    if(!keep_ts2){
      ts2[a:b] <- ts1  
    } 
    out <- ts2
  } else if(ts1s > ts2s & ts1e > ts2e){
    out <- c(ts2[1:(which(time(ts2) == ts1s)-1)],ts1)
    out <- ts(out,start = ts2s,frequency = freq)
  } else {
    stop('Case not covered, try switching ts1, ts2.')
  }
  
  
  out
  
  
}

