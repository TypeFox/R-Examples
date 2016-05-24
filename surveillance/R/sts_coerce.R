######################################################################
#Conversion between ts objects and sts objects
######################################################################

### Convert a simple "ts" object to an "sts" object

setAs(from = "ts", to = "sts", def = function (from) {
    ## Extract frequency and start from the "ts" object
    freq <- frequency(from)
    start <- start(from)
    if (length(start) == 1)
        stop("could not convert time series start() to (year, index) form")
    
    ## Remove "tsp" attribute and "ts"/"mts" class
    tsp(from) <- NULL
    ## "tsp<-"(x,NULL) is documented to also remove "ts" and "mts" classes
    ## but in R < 3.3.0, it did not remove "mts" (see PR#16769)
    from <- unclass(from)
    
    ## Create the sts object
    .sts(observed = from, start = start, freq = freq)
})


### Convert an "sts" object to a simple "ts" object

as.ts.sts <- function (x, ...)
{
    ts(data = x@observed, start = x@start, frequency = x@freq)
}
setAs(from = "sts", to = "ts", def = function (from) as.ts.sts(from))


### Convert an "sts" object to an eXtensible Time Series "xts"

as.xts.sts <- function (x, order.by = epoch(x, as.Date = TRUE), ...)
{
    if (!missing(order.by) || x@freq %in% c(52, 365)) {
        xts::xts(x = x@observed, order.by = order.by, ...)
    } else {
        ## frequencies 4 and 12 are nicely handled by the as.xts.ts method
        xts::as.xts(as.ts.sts(x), ...)
    }
}


######################################################################
#Method to convert sts object to a data frame suitable for regression
#Params:
# row.names - from generic R function
# optional  - from generic R function
# freqByWeek -- if TRUE use information in week (supposed to be Dates)
#               to freq (e.g. used for regression model)
######################################################################

setMethod("as.data.frame", signature(x="sts"), function(x,row.names = NULL, optional = FALSE, ...) {
  #Convert object to data frame and give names
  res <- data.frame("observed"=x@observed, "epoch"=x@epoch, "state"=x@state, "alarm"=x@alarm,"population"=x@populationFrac)

  if (ncol(x) > 1) {
    colnames(res) <-  c(paste("observed.",colnames(x@observed),sep=""),"epoch",
                        paste("state.",colnames(x@observed),sep=""),
                        paste("alarm.",colnames(x@observed),sep=""),
                        paste("population.",colnames(x@observed),sep=""))
  } else {
      colnames(res) <-  c("observed","epoch","state","alarm","population")
  }
  
  #Add a column denoting the number of week
  if (x@epochAsDate) {
    #Convert to date
    date <- as.Date(x@epoch, origin="1970-01-01")
    epochStr <- switch( as.character(x@freq), 
                       "12" = "%m",
                       "52" =  "%V",
                       "365" = "%j")
                       
    #Find out how many epochs there are each year
    years <- unique(as.numeric(formatDate(date,"%Y")))
    dummyDates <- as.Date(paste(rep(years,each=6),"-12-",26:31,sep=""))
    maxEpoch <- tapply( as.numeric(formatDate(dummyDates, epochStr)), rep(years,each=6), max)
    #Assign this to result
    res$freq <- maxEpoch[pmatch(formatDate(date,"%Y"),names(maxEpoch),duplicates.ok=TRUE)]
    res$epochInPeriod <- as.numeric(formatDate(date,epochStr)) / res$freq
  } else {
    #Otherwise just replicate the fixed frequency
    res$freq <- x@freq
    res$epochInPeriod <- x@epoch %% res$freq
  }
  
  return(res)
})


