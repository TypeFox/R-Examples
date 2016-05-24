Ping <- function(url, pingArgs='', warn=NA,
                 show.output.on.console=FALSE){
##
## 1.  get host
##
    urlSplit0 <- strsplit(url, '://')[[1]]
    urlS0 <- urlSplit0[min(2, length(urlSplit0))]
    URL <- strsplit(urlS0, '/')[[1]][1]
##
## 2.  construct ping command
##
    pingCmd <- paste('ping', pingArgs, URL)
##
## 3.  issue ping command
##
    if(is.na(warn)){
        rawResults <- system(pingCmd, intern=TRUE,
               show.output.on.console=show.output.on.console)
    } else {
        op <- options(warn=warn)
        rawResults <- system(pingCmd, intern=TRUE,
               show.output.on.console=show.output.on.console)
        options(op)
    }
##
## 4.  could not find host?
##
    cnf <- grep('Ping request could not find host', rawResults)
    if(length(cnf)>0){
        return(list(rawResults=rawResults,
                    rawNumbers=numeric(0),
                    counts=c(sent=0, received=0, lost=0),
                    p.lost=NA,
                    stats=c(min=NA, avg=NA, max=NA, mdev=NA) ) )
    }
##
## 5.  raw times
##
    time. <- grep('time=', rawResults, value=TRUE)
    if(length(time.)>0){
        timesel <- strsplit(time., 'time=')
        timechar <- sapply(timesel, '[', 2)
        timech. <- strsplit(timechar, 'ms')
        timech1 <- sapply(timech., '[', 1)
        rawNumbers <- as.numeric(timech1)
        min. <- min(rawNumbers)
        avg. <- mean(rawNumbers)
        max. <- max(rawNumbers)
        mdev <- sd(rawNumbers)
    } else{
        rawNumbers <- numeric(0)
        min. <- avg. <- max. <- mdev <- NA
    }
##
## 6.  timed out
##
    out. <- grep('timed out', rawResults, value=TRUE)
##
## 7.  counts
##
    rcvd <- length(rawNumbers)
    lost <- length(out.)
    counts <- c(sent=rcvd+lost, received=rcvd, lost=lost)
##
## 8.  p.lost
##
    p.lost <- lost/counts[1]
##
## 9.  stats
##
    stats <- c(min=min., avg=avg.,
               max=max., mdev=mdev )
##
## 10.  Done
##
    list(rawResults=rawResults, rawNumbers=rawNumbers,
         counts=counts, p.lost=p.lost, stats=stats)
}
