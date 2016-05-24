#' Average an Oklahoma Mesonet time series data frame
#' 
#' Summarize an Oklahoma Mesonet time series (MTS) data 
#' frame returned by \code{\link{okmts}}. Summary can be by hour, day, 
#' month, or year. Precipitation (variable RAIN) is returned as totals within 
#' the given time period.
#'
#' @param mts data frame returned by \code{okmts}.
#' @param by character string indicating time period to average over. May
#'  include "hour", "day", "month", or "year".
#' @param metric function to summarize with. Default is "mean" (average), but
#' may also include "min", "max", and "sd" for minimum, maximum, and 
#' standard deviation, respectively.

#' @export

#' @return A data frame summarizing Mesonet measurements by station and given
#' time period.

#' @examples
#' \dontrun{
#' ## Retrieve Bessie station MTS files for 00:00 Jun 01, 1997
#' ## through 23:55 Oct 31, 1997
#' bess.mts <- okmts(begintime="1997-06-01 00:00:00",
#'  endtime="1997-10-31 23:55", station="bess")
#'
#' ## Average MTS data by day.
#' bess.mts.avg  <- avgokmts(bess.mts, by="day")
#' }


avgokmts <- function(mts, by, metric="mean") {
  ## Averages MTS data frame by hour, day, month, or year
  ## Arguments:
  ##  mts: MTS data frame provided by okmts()
  ##  by: character, indicating time period to average over:
  ##    hour, day, month, and/or year
  ## Returns: data frame
  
  ## change by to lowercase
  by <- tolower(by)
  
  if(length(by)>1) {
    stop(paste("Only one grouping (by) variable allowed"))
  }
  
  ## check by for appropriate grouping
  if(any(by %in% c("hour", "day", "month", "year"))==FALSE) {
    stop(c("Grouping (by) variable must be hour, day, month, or year"))
  }
  
  ## check metric for mean, max, or min
  if(any(metric %in% c("mean","max","min", "sd"))==FALSE) {
    stop(c("metric must be mean, max, min, or sd"))
  }
  
  ## set list for grouping variables
  by.list <- vector(mode="list", length=length(by)+2)
  
  ## set first grouping to station, identified by mts$STID
  by.list[[1]] <- mts$STID
  ## set second grouping to station number, identified by mts$STNM
  by.list[[2]] <- mts$STNM
  
  ## set grouping variables
  if(by=="hour") {
    by.list[[3]] <- format(mts$TIME, "%H")
    by.list[[4]] <- format(mts$TIME, "%d")
    by.list[[5]] <- format(mts$TIME, "%m")
    by.list[[6]] <- format(mts$TIME, "%Y")
    names(by.list) <- c("STID", "STNM", "HOUR", "DAY", "MONTH", "YEAR")
  } else if (by=="day") {
      by.list[[3]] <- format(mts$TIME, "%d")
      by.list[[4]] <- format(mts$TIME, "%m")
      by.list[[5]] <- format(mts$TIME, "%Y")
      names(by.list) <- c("STID", "STNM", "DAY", "MONTH", "YEAR")
  } else if(by=="month") {
    by.list[[3]] <- format(mts$TIME, "%m")
    by.list[[4]] <- format(mts$TIME, "%Y")
    names(by.list) <- c("STID", "STNM","MONTH", "YEAR")
  } else if(by=="year") {
    by.list[[3]] <- format(mts$TIME, "%Y")
    names(by.list) <- c("STID", "STNM", "YEAR")
  }

  ## variables to average
  avg.var <- colnames(mts)[names(mts)!="STID" & names(mts)!="STNM" 
                          & names(mts)!="TIME" ]
  
  ## calculate averages based on grouping variables
  mts.avg <- aggregate(mts[,avg.var], by=by.list, FUN=metric, na.rm=T)
  
  ## if avg.var is only one variable, it is returned as "x" in mts.avg
  ## change "x" to appropriate name
  if(length(avg.var)==1) {
    colnames(mts.avg)[colnames(mts.avg)=="x"] <- avg.var
  }
  
  ## calculate rain average
  if(any(colnames(mts) %in% "RAIN")) {
    mts.avg$RAIN  <- totalprecip(mts, by)
  }
  return(mts.avg)
}