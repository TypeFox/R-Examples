#' Retrieve Oklahoma Mesonet climatological data
#' 
#' Retrieve \href{http://www.mesonet.org/}{Oklahoma Mesonet} 
#' time series (MTS) data for a given time period and station. Alternatively, 
#' if station is omitted and latitude and longitude are given, retrieve  
#' MTS data for the closest operating station during the given time period.
#'
#' The Oklahoma Mesonet is a network of 
#' automated climate monitoring stations throughout Oklahoma, USA. 
#' Data collection began January 01, 1994; as of November 2014, there are 120 
#' active stations. Measurements are recorded every five minutes and sent to a 
#' central facility for verification and quality control by the Oklahoma 
#' Climatological Survey.
#'
#' Data access may be restricted by organization and/or location. Please refer 
#' to and follow policies found within the
#' \href{http://www.mesonet.org/index.php/site/about/data_access_and_pricing}{Oklahoma 
#' Mesonet Data Access Policy}. The authors and maintainers of \code{okmesonet}
#' assume no responsibility for the use or misuse of \code{okmesonet}.
#'
#' The objects used to define the time period for \code{okmts} can be either
#' character strings or POSIXct objects. Character strings should be in the 
#' format "\code{2009-09-08 09:05}" or "\code{2005-12-13 00:00:00}". POSIXct 
#' objects need to have a time zone specified; \code{okmts} converts time 
#' zones appropriately to download correct MTS data.
#'
#' Four letter Mesonet station identifier can be found in 
#' \code{\link{okstations}} or on the 
#' \href{http://www.mesonet.org/}{Oklahoma Mesonet} website.
#'
#' Available Mesonet variables and units are described in the 
#' \href{http://www.mesonet.org/index.php/site/about/mdf_mts_files}{MDF/MTS 
#' Files} webpage, 'Parameter Description' 
#' \href{http://www.mesonet.org/files/parameter_description_readme.pdf}{readme}
#' file, or \href{http://www.mesonet.org/wiki/Public:MDF_Format}{MTS 
#' specification}. Multiple variables can be retrieved by combining values into
#' a vector, e.g. \code{c("TAIR", "RELH")}. \code{"ALL"} indicates all 
#' available variables.
#'
#' Time records of Oklahoma MTS data are stored in Coordinated Universal Time
#' (UTC or GMT). To easily convert to local Oklahoma time, \code{localtime=TRUE}
#' indicates that times used to define the time period (\code{begintime} and
#' \code{endtime}) are local Oklahoma time.
#' Time zone conversion is done internally, and accounts for Daylight Savings
#' Time (as reliably as R can; see \link{timezone}).
#' \code{localtime=TRUE} will also direct \code{okmts} to output in local 
#' Oklahoma time. \code{localtime=FALSE} indicates that UTC is used for 
#' both \code{begintime} and \code{endtime}; output is also UTC. If time 
#' inputs are of POSIXct class, \code{localtime} only affects time output.
#'
#' Missing values are stored as negative integer codes and can be converted to
#' NA with the \code{missingNA} parameter. Missing value descriptions can be 
#' found in the
#' \href{http://www.mesonet.org/index.php/site/about/mdf_mts_files}{MDF/MTS 
#' Files} webpage.
#' 
#' The use of multiple cores can speed up data retrieval for lengthy time 
#' periods. \code{mcores} specifies the number of cores to be used. 
#' \code{mcores=TRUE} will direct \code{okmts} to use the 
#' number cores less one in the current machine (determined by 
#' \code{\link[parallel]{detectCores}-1}).
#'
#' To prevent repeated retrieval of frequently used data, the data frame
#' returned by \code{okmts} can be saved (e.g. \code{\link{save}}) or 
#' written to a file (e.g. \code{\link{write.table}}).
#'
#' @param begintime character string or POSIXct object. Start time of time 
#' period. Character strings must be formatted as 'YYYY-MM-DD HH:MM:SS'.
#' @param endtime character string or POSIXct object. End time of time 
#' period. Character strings must be formatted as 'YYYY-MM-DD HH:MM:SS'.
#' @param station character string. Four letter Mesonet station identifier. 
#' See 'Details'.
#' @param lat numeric. latitude of point of interest in decimal degrees.
#' @param lon numeric. longitude of point of interest in decimal degrees.
#' @param variables character string. Mesonet variables to retrieve. See 'Details'.
#' @param localtime logical; if \code{TRUE}, input and output time is local to
#'  Oklahoma. If \code{FALSE}, input and output time is Coordinated Universal 
#'  Time (UTC or GMT). See 'Details'.
#' @param missingNA logical; if \code{TRUE}, missing values are replaced with
#'  NA. See 'Details'.
#' @param mcores integer or logical; use \emph{n} cores for file retrieval.
#'  See 'Details'.

#' @export

#' @seealso \code{\link{avgokmts}} to summarize MTS data.

#' @return A data frame with values from MTS files for the given station, time 
#' period, and desired variables. Time values for each measurement are returned 
#' as POSIXct class; time zone is determined by \code{localtime}.

#' @examples
#' \dontrun{
#' ## Retrieve Bessie station MTS files for 00:00 Jun 01, 1997
#' ## through 23:55 Oct 31, 1997
#' bess.mts <- okmts(begintime="1997-06-01 00:00:00",
#'  endtime="1997-10-31 23:55", station="bess")
#'
#' ## Use POSIXct class to retrieve Medicine Park station air
#' ## temperature for 09:30 through 20:30 Aug 12, 2004
#' ## Set times, using 'America/Chicago' for Oklahoma time zone
#' medi.time <- c(as.POSIXct("2004-08-12 09:30", tz="America/Chicago"),
#'  as.POSIXct("2004-08-12 20:30", tz="America/Chicago"))
#' medi.mts <- okmts(begintime=medi.time[1], endtime=medi.time[2],
#'  station="medi", variables=c("TAIR", "RELH"))
#'
#' ## Download all data for 2001 for station closest to 
#' ## 36.575284 latitude, -99.478455 longitude, using multiple cores
#' stn.mts <- okmts(begintime="2001-01-01 00:00:00", 
#'  endtime="2001-12-31 23:55:00", lat=36.575284, lon=-99.478455, mcores=T)
#'
#' ## Retrieve Idabel station MTS data for 00:00 through 12:00 UTC (GMT)
#' ## Nov 23, 2003
#' ## Time values are returned in UTC
#' idab.mts <- okmts(begintime="2003-11-23 00:00:00", 
#'  endtime="2003-11-23 12:00:00", station="idab", localtime=F)
#'
#' ## Combine air temperature with bison movement data.
#' ## Retrieve Foraker station MTS files for 00:00 Jan 31, 2011 
#' ## through 15:00 Feb 05, 2011
#' fora.mts <- okmts(begintime="2011-01-31 00:00:00", 
#'  endtime="2011-02-05 15:00:00", station="fora")
#' ## Round bison timestamp down to five minute mark
#' bison$newtime <- round(bison$timestamp, "min")
#' bison$newtime$min <- as.integer(format(bison$newtime, "%M")) %/% 5 * 5
#' bison$newtime <- as.POSIXct(bison$newtime)
#' ## Add Foraker station air temperature to bison data
#' bison$TAIR <- fora.mts$TAIR[match(bison$newtime, fora.mts$TIME)]
#' }

okmts <- function(begintime, endtime, station=NULL, lat=NULL, lon=NULL, 
                  variables="ALL", localtime=TRUE, missingNA=TRUE,
                  mcores=FALSE) {
  ## Gets Mesonet MTS file from Mesonet homepage
  ## Arguments:
  ##  begintime: beginning date,given as 'YYYY-MM-DD HH:MM:SS'
  ##  endtime: end date, given as 'YYYY-MM-DD HH:MM:SS'
  ##  station: four letter character ID for Mesonet station
  ##  lat: latitude of point location, decimal degrees
  ##  lon: longitude of point locaiton, decimal degrees
  ##  variables: variables to retrieve
  ##  localtime: logical to indicate the use of Oklahoma local time, else 
  ##    use GMT
  ##  missingNA: logical to indicate replace missing values with NA
  ##  mcores: logical to indicate use of foreach and multiple cores
  ## Returns: dataframe of MTS files
  
  ## check to see if station information is available
  if(exists("okstations")==F || nrow(okstations)<50) {
    path.geoinfo <- paste("http://www.mesonet.org/index.php/api/siteinfo/",
                          "from_all_active_with_geo_fields/format/csv/", 
                          sep ="")
    stop.msg <- paste("Oklahoma Mesonet station list unavailable or", 
                      "incomplete. Check", path.geoinfo,
                      "for connectivity and run", 
                      sQuote("okstations  <- updatestn()"), 
                      "to update station list")
    stop(stop.msg)
  }
  
  ## check to see if begintime and endtime are of class character or POSIXct 
  ## set *.local and *.gmt appropriately
  if(is.character(begintime)==T & is.character(endtime)==T) {
    if(localtime==T) {
      ## convert character timestamp to POSIXct, timezone America/Chicago
      begintime.local <- as.POSIXct(begintime, tz="America/Chicago")
      endtime.local  <- as.POSIXct(endtime, tz="America/Chicago")
      ## Convert to timezone GMT for file retrieval
      begintime.gmt <- as.POSIXct(format(begintime.local, tz="GMT"), tz="GMT")
      endtime.gmt <- as.POSIXct(format(endtime.local, tz="GMT"), tz="GMT")
    } else {
      ## convert character timestamp to POSIXct, timezone GMT for file retrieval
      begintime.gmt <- as.POSIXct(begintime, tz="GMT")
      endtime.gmt  <- as.POSIXct(endtime, tz="GMT")
      ## convert begintime.gmt to timezone America/Chicago for continuity
      begintime.local <- as.POSIXct(format(begintime.gmt, tz="America/Chicago"),
                                    tz="America/Chicago")
      endtime.local <- as.POSIXct(format(endtime.gmt, tz="America/Chicago"),
                                  tz="America/Chicago")
    }
  } else if(any(class(begintime)=="POSIXct") & any(class(endtime)=="POSIXct")) {
    if(localtime==T) {
      ## convert POSIXct timestamp to timezone America/Chicago
      ## used for subsetting desired variables
      begintime.local <- as.POSIXct(format(begintime, tz="America/Chicago"),
                                    tz="America/Chicago")
      endtime.local <- as.POSIXct(format(endtime, tz="America/Chicago"),
                                  tz="America/Chicago")
      ## convert POSIXct timestamp to timezone GMT for file retrieval
      begintime.gmt <- as.POSIXct(format(begintime, tz="GMT"), tz="GMT")
      endtime.gmt <- as.POSIXct(format(endtime, tz="GMT"), tz="GMT")
    } else {
      ## convert POSIXct timestamp to timezone GMT for file retrieval
      begintime.gmt <- as.POSIXct(format(begintime, tz="GMT"), tz="GMT")
      endtime.gmt <- as.POSIXct(format(endtime, tz="GMT"), tz="GMT")
      ## convert begintime.gmt to timezone America/Chicago for continuity
      begintime.local <- as.POSIXct(format(begintime.gmt, tz="America/Chicago"),
                                    tz="America/Chicago")
      endtime.local <- as.POSIXct(format(endtime.gmt, tz="America/Chicago"),
                                  tz="America/Chicago")
    }
  } else {
    ## if not character or POSIXct, stop and give error message
    stop(paste("begintime and endtime must both be entered as",
               dQuote("YYYY-MM-DD HH:MM:SS"), "or a POSIXct class."))
  }
  
  ## if station is NULL and lat and long are given, retrieve closest station
  ## with nearstn()
  if(is.null(station)==T & is.numeric(lat)==T & is.numeric(lon)==T) {
    ## check lat/long coordinates first to make sure not too far away
    if(lat>38 | lon<c(-104) | lat<33 | lon>c(-94))
      stop(paste(lat, "Latitude,", lon, "Longitude is too far away from",
                 "Oklahoma for meaningful data."))
    station <- nearstn(pnt.lon=lon, pnt.lat=lat, startdate=begintime.local, 
                       enddate=endtime.local)
  }
  
  ## check to see if station is a true station
  if(toupper(station) %in% okstations$Identifier==F) {
    stop.msg <- paste("Station identifier is incorrect.",
                      "Please check", sQuote("okstations"), 
                      "or http://www.mesonet.org/ for correct four letter", 
                      "identifier.")
    stop(stop.msg)
  }
  
  ## verify begintime is before endtime
  if(begintime.gmt > endtime.gmt) {
    stop.msg <- paste("Parameter begintime (", begintime, ") is after ",
                      "endtime (", endtime, ").", sep="")
    stop(stop.msg)
  }
  
  ## check to see if begintime and endtime are before station commission date
  comm.date.local <- okstations$Commissioned[match(toupper(station),
                                             okstations$Identifier)]
  comm.date.gmt <- as.POSIXct(format(comm.date.local, tz="GMT"), tz="GMT")
  if(begintime.gmt < comm.date.gmt | endtime.gmt < comm.date.gmt){
    stop.msg <- paste("Parameters begintime or endtime are before station was", 
                      " commissioned (", 
                      format(comm.date.local, "%Y-%m-%d"), 
                      ").", sep="")
    stop(stop.msg)
  }
  
  ## check to see if begintime and endtime are before station decommissioned 
  ## date
  decomm.date.local <- okstations$Decommissioned[match(toupper(station),
                                                   okstations$Identifier)]
  decomm.date.gmt <- as.POSIXct(format(decomm.date.local, tz="GMT"), tz="GMT")
  if(begintime.gmt > decomm.date.gmt | endtime.gmt > decomm.date.gmt){
    stop.msg <- paste("Parameters begintime or endtime are after station was", 
                      " decommissioned (", 
                      format(decomm.date.local, "%Y-%m-%d"), 
                      ").", sep="")
    stop(stop.msg)
  }
  
  ## available Mesonet variables
  mtsvariables <- c("STID", "STNM", "RELH", "TAIR", "WSPD", "WVEC", "WDIR", "WDSD", 
                 "WSSD", "WMAX", "RAIN", "PRES", "SRAD", "TA9M", "WS2M", "TS10", 
                 "TB10", "TS05", "TB05", "TS30", "TR05", "TR25", "TR60", "TR75",
                 "ALL")
  
  ## convert variables to uppercase
  variables <- toupper(variables)
    
  ## check to see if variables matches available variables
  if(all(variables %in% mtsvariables)==FALSE) {
    stop.msg <- paste("Desired variables do not match available variables. See",
                      "http://www.mesonet.org/index.php/site/about/mdf_mts_files",
                      "for available variables.")
    stop(stop.msg) 
  }
  
  ## if variables contains "ALL", remove anything else
  if(any(variables %in% "ALL")==TRUE) variables <- "ALL"
  
  ## convert station to lowercase
  station <- tolower(station)
  
  ##  sequence GMT days from begin to end for file retrieval
  dates.gmt <- seq.POSIXt(trunc(begintime.gmt, units="days"),
                          trunc(endtime.gmt, units="days"), by="days")
  
  ## verify access to MTS files
  if(length(dates.gmt)<10) nverify <- length(dates.gmt) else nverify <- 10
  verifyMTS <- sapply(dates.gmt[1:nverify], FUN=verifymts, station=station)
  if(all(verifyMTS==F)) {
    stop.msg <- paste("Access to the first ", nverify, " MTS files is ",
                      "unavailable. Please check ",
                      "http://www.mesonet.org/index.php/dataMdfMts/",
                      "dataController/getFile/", 
                      format.POSIXct(dates.gmt[1], format="%Y%m%d"), station, 
                      "/mts/TEXT/ ", "for access.", sep="")
    stop(stop.msg)
  }
  
  ## create empty lists
	all.MTS <- vector(mode="list", length=length(dates.gmt))
  
  ## use multiple cores if indicated by mcores=T
  if(mcores==T | is.numeric(mcores)==T) {
    if(mcores==T) ncores <- detectCores()-1
    if(is.numeric(mcores)==T) ncores=round(mcores)
    if(.Platform$OS.type=="unix") {
      all.MTS <- mclapply(dates.gmt, FUN=retrievemts, station=station,
                          variables=variables, missingNA=missingNA,
                          mc.cores=ncores)
    } else if(.Platform$OS.type=="windows") {
      c1 <- makeCluster(getOption("cl.cores", ncores))
      all.MTS <- parLapply(c1, dates.gmt, fun=retrievemts, station=station, 
                           variables=variables, missingNA=missingNA)
      stopCluster(c1)
    } } else {
    all.MTS <- lapply(dates.gmt, FUN=retrievemts, station=station, 
                      variables=variables, missingNA=missingNA)
  }
  
  ##  If localtime==T, convert back to CST/CDT and subset to begin/end time
  if(localtime==T) {
    all.MTS <- lapply(all.MTS,
                      function(x) {
                        x$TIME <- 
                          as.POSIXct(format(x$TIME, tz="America/Chicago"),
                                     tz="America/Chicago")
                        return(x)
                      })
    ##  Subset data according to begin and end time
    all.MTS[[1]] <- subset(all.MTS[[1]], TIME>=begintime.local & 
                           TIME<=endtime.local)
    all.MTS[[length(all.MTS)]] <- subset(all.MTS[[length(all.MTS)]],
                                         TIME>=begintime.local & 
                                           TIME<=endtime.local)
  } else {
    ##  Subset data according to begin and end time
    all.MTS[[1]] <- subset(all.MTS[[1]], TIME>=begintime.gmt & 
                           TIME<=endtime.gmt)
    all.MTS[[length(all.MTS)]] <- subset(all.MTS[[length(all.MTS)]],
                                         TIME>=begintime.gmt & 
                                           TIME<=endtime.gmt)
  }
  
  ## use ldply from plyr package to return list as dataframe
  return(ldply(all.MTS))
}
