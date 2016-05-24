#' Get raw data from Markit website.
#' 
#' \code{get_raw_markit} downloads the rates zip file from Markit website,
#' unzips and parses the XML
#' 
#' @param date Date type, is the CDS pricing date.
#' @param currency numeric, is the currency that the CDS is traded in.
#'   
#' @importFrom utils unzip
#'   
#' @return a data frame that contains CDS pricing date, currency, interest
#'  rate expiry and interest rate. The data frame is created with data
#'  from Markit website

get_raw_markit <- function(date, currency){
  
  stopifnot(inherits(date, "Date") & is.character(currency))
  stopifnot(currency %in% c( "USD", "EUR", "JPY"))
  
  ## CDS for Trade Date will use rates from Trade Date - 1 
  
  date <- date - 1
  
  ## 0 is Sunday, 6 is Saturday
  
  dateWday <- as.POSIXlt(date)$wday
  
  ## change date to the most recent weekday if necessary i.e. change date to 
  ## Friday if the day we are pricing CDSs is Saturday or Sunday
  
  if (dateWday == 0){
    date <- date - 2
  } else if (dateWday == 6) {
    date <- date - 1
  }
  
  ## convert date to numeric
  
  dateInt <- format(date, "%Y%m%d")
  
  ## markit rates URL
  
  ratesURL <- paste("https://www.markit.com/news/InterestRates_",
                    currency, "_", dateInt, ".zip", sep ="")
  
  ## define an internal function here for URL handling
  
  rates.URL <- function(URL, verbose = FALSE){ 
    
    tf <- tempfile()
    td <- tempdir()
    
    f <- CFILE(tf, mode="wb")
    a <- tryCatch(curlPerform(url = URL,
                              writedata = f@ref, noprogress=TRUE,
                              verbose = FALSE,
                              ssl.verifypeer = FALSE),
                  error = function(e) {
                    return("Rates data not available at markit.com")
                  })
    close(f)
    if (class(a) == "character"){
      return(a)
    } else {
      
      ## add suppressWarnings here because otherwise the stop message of Internet
      ## Connection Problem looks confusing.
      
      files <- suppressWarnings(utils::unzip(tf , exdir = td))
      
      ## the 2nd file of the unzipped directory contains the rates info
      
      ## add clearer message of Internet Connection Problem
      
      if(length(files) == 0){
        stop(
          "Either an internet connection problem occurs or the data
        of the required date are unavailable on the internet. Try 
        again later.")
      }
      
      doc <- xmlTreeParse(files[grep(".xml", files)], getDTD = F)
      return(xmlRoot(doc))
    }
  }
  
  
  ## XML file from the internet, which contains the rates data
  
  xmlParsedIn <- rates.URL(ratesURL)
  
  ## rates data extracted from XML file
  
  rates <- xmlSApply(xmlParsedIn, function(x) xmlSApply(x, xmlValue))
  
  ## extracts the 'M' or 'Y' of the expiry and stores it in curveRates
  
  curveRates <- c(rates$deposits[names(rates$deposits) == "curvepoint"],
                  rates$swaps[names(rates$swaps) == "curvepoint"])
  
  ## split the numbers from the 'M' and 'Y'
  
  x <- do.call(rbind, strsplit(curveRates, split = "[MY]", perl = TRUE))
  rownames(x) <- NULL
  x <- cbind(x, "Y")
  
  ## attacg M to money month rates
  
  x[1: (max(which(x[,1] == 1)) - 1), 3] <- "M"
  
  ## data frame with Interest Rates, maturity, type, expiry
  
  ratesx <- data.frame(expiry = paste(x[,1], x[,3], sep = ""),
                       matureDate = substring(x[,2], 0, 10),
                       rate = substring(x[,2], 11),
                       
                       ## if maturity is 1Y, it is of type M
                       
                       type = c(rep("M", sum(names(rates$deposits) == "curvepoint")),
                                rep("S", sum(names(rates$swaps) == "curvepoint"))))
  
  ratesx$expiry <- as.character(ratesx$expiry)
  
  return(ratesx)
}
