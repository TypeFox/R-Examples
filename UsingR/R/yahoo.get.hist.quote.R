#######################################################
### adapted from tseries, RMetrics
### author: Daniel Herlemont <dherlemont@yats.com>
#######################################################

yahoo.get.hist.quote <-
  function(
           instrument = "^gspc"
           ,destfile = paste(instrument,".csv",sep="")
           ,start, end, quote = c("Open", "High", "Low", "Close"),
           adjusted=TRUE, download=TRUE,
           origin = "1970-01-01", compression = "d")  {
    
#######################################################
### read local quotes file
### adapted from tseries, RMetrics
### author: Daniel Herlemont <dherlemont@yats.com>
#######################################################


  yahoo.get.hist.quote.local <-
    function (instrument = "^gspc"
              ,destfile = paste(instrument,".csv",sep="")
              ,start, end, quote = c("Open", "High", "Low", "Close"),
              adjusted=TRUE, origin = "1970-01-01") {

      if(missing(start)) start <- "1991-01-02"
      if(missing(end))
        end <- format(Sys.time() - 86400, "%Y-%m-%d")
      start <- as.POSIXct(start, tz = "GMT")
      end <- as.POSIXct(end, tz = "GMT")
                                        # !!! patch, remove all download stuf
      nlines <- length(count.fields(destfile, sep = "\n"))
      if(nlines == 1) {
        stop(paste("No data available for", instrument))
      }
      
      x <- read.table(destfile, header = TRUE,
                      sep = ",", as.is = TRUE, fill = TRUE)
      x <- na.omit(x)
      
                                        # !!! do adjustements
      if (adjusted) {
        adjust=x[,7]/x[,5] # factor
        x[,2]=x[,2]*adjust # 
        x[,3]=x[,3]*adjust
        x[,4]=x[,4]*adjust
        x[,5]=x[,7]					# close = adjusted close
        x[,6]==x[,6]/adjust # divide for volume
      }
    
      ##  remove file
      unlink(destfile)
    
      nser <- pmatch(quote, names(x)[-1]) + 1
      if(any(is.na(nser)))
        stop("This quote is not available")
      n <- nrow(x)
    
      ## Yahoo currently formats dates as '26-Jun-01', hence need C
      ## LC_TIME locale for getting the month right.
      lct <- Sys.getlocale("LC_TIME")
      Sys.setlocale("LC_TIME", "C")
      on.exit(Sys.setlocale("LC_TIME", lct))
      
      dat <- gsub(" ", "0", as.character(x[, 1])) # Need the gsub?
      dat <- as.POSIXct(strptime(dat, "%d-%b-%y"), tz = "GMT")
      ## quiet this function
      ##      if(dat[n] != start)
      ##        cat(format(dat[n], "time series starts %Y-%m-%d\n"))
      ##      if(dat[1] != end)
      ##        cat(format(dat[1], "time series ends   %Y-%m-%d\n"))
      jdat <-
        unclass(julian(dat, origin = as.POSIXct(origin, tz = "GMT")))
      ## We need unclass() because 1.7.0 does not allow adding a number
      ## to a "difftime" object.
      ind <- jdat - jdat[n] + 1
      y <- matrix(NA, nrow = max(ind), ncol = length(nser))
      y[ind, ] <- as.matrix(x[, nser, drop = FALSE])
      colnames(y) <- names(x)[nser]
      return(ts(y, start = jdat[n], end = jdat[1]))
    }
  
  
  
#######################################################
### download quotes file
### adapted from tseries, RMetrics
### author: Daniel Herlemont <dherlemont@yats.com>
#######################################################
  yahoo.get.hist.quote.download=function (
    instrument = "^gspc"
    ,destfile = paste(instrument,".csv",sep="")
    ,start, end, origin = "1970-01-01", compression = "d") 
    {
      if (missing(start)) 
        start <- "1991-01-02"
      if (missing(end)) 
        end <- format(Sys.time() - 86400, "%Y-%m-%d")
      start <- as.POSIXct(start, tz = "GMT")
      end <- as.POSIXct(end, tz = "GMT")
      url <- paste("http://chart.yahoo.com/table.csv?s=", instrument, 
                   format(start,
                          paste(
                                "&a=", as.character(as.numeric(format(start,"%m"))
                                                    - 1),
                                "&b=%d&c=%Y", sep = "")),
                   format(end,
                          paste("&d=",
                                as.character(as.numeric(format(end, "%m")) - 1), "&e=%d&f=%Y", sep = "")),
                   "&g=", 
                   compression,
                   "&q=q&y=0&z=", instrument,
                   "&x=.csv", 
                   sep = "")
      
      status <- download.file(url, destfile)
      if (status != 0) {
        unlink(destfile)
        stop(paste("download error, status", status))
      }
    }
  
  
  if (download) {
    yahoo.get.hist.quote.download(instrument=instrument,destfile=destfile
                                  ,start=start, end=end, origin = origin,
                                  compression=compression) 
  }
  yahoo.get.hist.quote.local(instrument=instrument,destfile=destfile
                             ,start=start, end=end
                             ,quote=quote,adjusted=adjusted
                             ,origin = origin) 
}

