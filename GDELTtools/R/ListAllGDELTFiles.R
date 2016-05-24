ListAllGDELTFiles <- function(as.of=Sys.Date() - 1) {
  compressed.files <- c(paste(1979:2005, ".zip", sep=""),
                        paste(format(seq(from=as.Date("2006-01-01"), to=as.Date("2013-03-01"), by="month"), "%Y%m"), ".zip", sep=""),
                        paste(format(seq(from=as.Date("2013-04-01"), to=as.Date(as.of), by="day"), "%Y%m%d"), ".export.CSV.zip", sep=""))
  data.files <- c(paste(1979:2005, ".csv", sep=""),
                  paste(format(seq(from=as.Date("2006-01-01"), to=as.Date("2013-03-01"), by="month"), "%Y%m"), ".csv", sep=""),
                  paste(format(seq(from=as.Date("2013-04-01"), to=as.Date(as.of), by="day"), "%Y%m%d"), ".export.CSV", sep=""))
  
  # Remove dates GDELT just doesn't have
  compressed.files <- setdiff(compressed.files, c("20140123.export.CSV.zip","20140124.export.CSV.zip","20140125.export.CSV.zip"))
  data.files <- setdiff(data.files, c("20140123.export.CSV","20140124.export.CSV","20140125.export.CSV"))
  
  
  return(list(compressed=compressed.files, data=data.files))
}