## Segment of code that returns list of all necessary GDELT files to cover the time span specified

# input startdate and enddate in the form "yyyy-mm-dd"

FileListFromDates <- function(start.date, end.date=start.date){  
  
  start.date <- as.Date(start.date)
  end.date <- as.Date(end.date)
  
  if(end.date < start.date) stop("end.date cannot be before start.date")
  if(start.date < as.Date("1979-01-01")) stop("start.date cannot be before 1979")
  if(end.date > Sys.Date()) stop("end.date cannot be in the future")
  
  out <- character(0)
  
  ## PART A
  ## make the yyyy list for 1979 through 2005
  
  if(start.date < "2006-01-01"){
    start.year <- as.numeric(format(start.date, "%Y"))
    end.year <- min(c(2005,as.numeric(format(end.date, "%Y"))))
    out <- c(out, paste(start.year:end.year, ".zip", sep=""))
  }
  
  # PART B
  # Make the yyyymm list for 2006 through March 2013
  
  if(start.date < as.Date("2013-04-01") & end.date >= as.Date("2006-01-01")) {
    
    if( start.date < as.Date("2006-01-01") ) {
      b.start.date <- as.Date("2006-01-01")
    } else {
      b.start.date <- as.Date(paste(format(start.date, "%Y-%m-"), "01", sep=""))
    }
    
    if( end.date < as.Date("2013-04-01") ) {
      b.end.date <- end.date
    } else {
      b.end.date <- as.Date("2013-03-31")
    }
    out <- c(out, paste(format(seq(from=b.start.date, to=b.end.date, by="month"), "%Y%m"), 
                        ".zip", sep=""))
  }
  
  # PART C
  # Make list of yyyymmdd starting April 1, 2013
  if( end.date >= as.Date("2013-04-01") ){
    if( start.date < as.Date("2013-04-01") ) {
      c.start.date <- as.Date("2013-04-01")
    } else {
      c.start.date <- start.date
    }
    out <- c(out, paste(format(seq(from=c.start.date, to=end.date, by="day"), "%Y%m%d"), 
                        ".export.CSV.zip", sep=""))
  }
  
  # Remove dates GDELT just doesn't have
  out <- setdiff(out, c("20140123.export.CSV.zip","20140124.export.CSV.zip","20140125.export.CSV.zip"))
  
  return( out )
}
