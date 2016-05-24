#' @title Download US Treasury yield curve data.
#' @description Downloads US Treasury yield curve data from the US Treasury web site.
#' @details Forms a query to submit for US Treasury yield curve data, posting this query to the US Treasury web site's data feed service.  By default the download includes data yield data for 12 products from January 1, 1990, some of which are NA during this span.  The caller can pass parameters to limit the query to a certain year or year and month, but the full download is not especially large.  The download data from the service is in XML format.  This function transforms that data into a numeric data frame with treasury product items (constant maturity yields for 12 kinds of bills, notes, and bonds) as columns and dates as row names. The function returns a list which includes an item for this data frame as well as query-related values for reference and the update date from the service.  The data frame can be used as-is or converted easily to a time series format such as \code{xts}. 
#' @param year the desired year number or NULL for all years (default)
#' @param month the desired month number or NULL for all months (default)
#' @param base the base URL for the data service, defaulting to \url{http://data.treasury.gov/feed.svc/DailyTreasuryYieldCurveRateData}.  If the month or year arguments are not NULL, then the function modifies this URL to parameterize the download request.
#' @param allowParallel whether to allow \code{ldply} to use a registered parallel cluster.  FALSE by default. 
#' @return Class type \code{ustyc} containing update date \code{updated}, dataframe \code{df}, \code{month}, \code{year}, and \code{query} elements.  The \code{query} element value is the string used to call the data service for download.
#' @export
#' @importFrom XML xmlParse xmlToList
#' @importFrom plyr ldply
#' @references \url{http://data.treasury.gov/feed.svc/DailyTreasuryYieldCurveRateData}
#' @seealso \url{http://cran.r-project.org/web/packages/FRBData/} for different interest rates and source. 
#' @examples
#' \dontrun{ 
#' xlist = getYieldCurve()
#' summary(xlist)
#' }
getYieldCurve <- function(year=NULL,
                          month=NULL,
                          base="http://data.treasury.gov/feed.svc/DailyTreasuryYieldCurveRateData",
                          allowParallel=FALSE
                          ) {
  #require(XML)
  #require(plyr)
  
  location <- base
  yloc <- mloc <- doc <- NULL
  yloc <- if(is.null(year)==FALSE) paste("year(NEW_DATE)%20eq%20",year,sep='')
  mloc <- if(is.null(month)==FALSE) paste("month(NEW_DATE)%20eq%20",month,sep='')
  
  # determine whether caller wants subset of data
  parameters <- ""
  if (is.null(yloc)==FALSE && is.null(mloc)==FALSE) {
    parameters = paste("?$filter=",mloc,"%20and%20",yloc,sep='')
  } else {
    if (is.null(yloc)==FALSE)
      parameters = paste("?$filter=",yloc,sep='')
    if (is.null(mloc)==FALSE)
      parameters = paste("?$filter=",mloc,sep='')
  }
  
  doc <- xmlParse(paste(location,parameters,sep=''))
  if (is.null(doc)) {
    warning(paste("Could not parse the location",location))
    return(NULL)
  }
  
  message("Download and parse complete.  Converting to list...")
  
  x <- xmlToList(doc)
  message("List conversion complete.  Converting to frame...")

  # save the updated time
  updated = x[[3]]
  
  # truncate first four elements and the last element
  x[1:4] <- NULL
  x[[length(x)]] <- NULL
  
  # field extraction function
  cy <- function(t,p) {
    if ("text" %in% names(p[[t]])) 
      p[[t]]$text
    else 
      NA
  }
  
  # list manipulator to produce a data frame
  y <- ldply(x, function(e) { 
    p <- e$content$properties
    q = sapply(names(p),cy,p)
    },
    .id="NEW_DATE",
    .parallel=allowParallel)

  # strip hours and sort by date
  y$NEW_DATE <- substring(y$NEW_DATE,1,10)
  y <- y[with(y,order(NEW_DATE)),]
  dates <- y$NEW_DATE
  
  # trim the columns, convert remainder to double, assign row names
  y <- data.frame(apply(y[,3:14],2,function(x) as.double(x)))
  rownames(y) <- dates
  
  message("Frame conversion complete.")

  # return a list with data frame and some other useful tags from fetch
  rv <- list(updated=updated,df=y,month=month,year=year,query=location)
  class(rv) <- "ustyc"
  rv
}

