read.testURLs <- function(file.='testURLresults.csv', ...){
##
## 1.  read(file., ...)
##
    dat <- read.csv(file., ...)
##
## 2.  parse Time
##
    if('Time' %in% names(dat)){
        tm <- dat$Time
        if(is.numeric(tm)){
#            if(require(fda)){
                Tm <- as.POSIXct1970(tm)
#            }
        } else {
            if(is.factor(tm))tm <- as.character(tm)
            if(is.character(tm)){
                Tm <- strptime(tm, '%a %b %d %H:%M:%S %Y', tz='GMT')
            } else TM <- tm
        }
        dat$Time <- Tm
    }
##
## 3.  Done
##
    class(dat) <- c('testURLs', 'data.frame')
    dat
}
