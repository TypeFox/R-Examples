#' Retrieve Oklahoma Mesonet current observations
#' 
#' Retrieve current observations from all
#' \href{http://www.mesonet.org/}{Oklahoma Mesonet} stations.
#'
#' See 
#' \href{http://www.mesonet.org/index.php/site/about/current_observations_csv}{Current Observations}
#' for variable and unit description.
#' Current observations are stored in 'American' units. \code{SI=TRUE} converts
#' to the 'International System of Units' (though Celsius is used instead of
#' Kelvin for temperature).
#'
#' A TIME variable of class POSIXct is added. \code{localtime=TRUE} returns 
#' TIME as local Oklahoma time; \code{localtime=FALSE} returns TIME as UTC.
#'
#' @param SI logical; if \code{TRUE}, convert observations to SI units. See 
#'  'Details'.
#' @param localtime logical; if \code{TRUE} output time is local to Oklahoma.
#'  If \code{FALSE}, output time is Coordinated Universal Time (UTC).

#' @export

#' @seealso \code{\link{okmts}}

#' @return A data frame with current observations.

#' @examples
#' \dontrun{
#' ## Retrieve current observations.
#' curobs <- okcurobs()
#' }


okcurobs <- function(SI=TRUE, localtime=TRUE) {
  ## Retreive current observations for all stations from Mesonet website
  ## Arguments:
  ##  SI: logical to indicate conversion to SI units
  ##  localtime: logical to indicate local Oklahoma time output
  ## Returns: dataframe of current observations, TIME variable added as POSIXct
  ##  class
  
  ## set path and read current observation file
  path <- paste("http://www.mesonet.org/data/public/mesonet/current/",
                "current.csv.txt", sep="")
  curobs <- read.csv(file=path,
                     colClasses=c("character", "character", "character",
                                  "numeric", "numeric", "integer", "integer",
                                  "integer", "integer", "integer",
                                  "integer", "integer", "integer", "integer",
                                  "integer", "character", "integer", "integer",
                                  "numeric", "integer", "integer", "numeric"))
  
  ## stitch timestamp together
  time.stitched <- paste(curobs$YR, sprintf("%02d", curobs$MO),
                         sprintf("%02d", curobs$DA),
                         sprintf("%02d", curobs$HR),
                         sprintf("%02d", curobs$MI), "-0600", sep="")
  ## add TIME variable
  curobs$TIME <- as.POSIXct(time.stitched, format="%Y%m%d%H%M%z")
  
  ## local Oklahoma time
  if (localtime==TRUE) {
    curobs$TIME <- as.POSIXct(format(curobs$TIME, tz="America/Chicago"),
                            tz="America/Chicago")
  } else {
    curobs$TIME <- as.POSIXct(format(curobs$TIME, tz="GMT"), tz="GMT")
  }
  
  ## if SI=TRUE, convert to SI units
  if(SI==T) {
    ## temperatures
    curobs$TAIR <- round((5/9) * (curobs$TAIR - 32), 1)
    curobs$TDEW <- round((5/9) * (curobs$TDEW - 32), 1)
    curobs$CHIL <- round((5/9) * (curobs$CHIL - 32), 1)
    curobs$HEAT <- round((5/9) * (curobs$HEAT - 32), 1)
    curobs$TMAX <- round((5/9) * (curobs$TMAX - 32), 1)
    curobs$TMIN <- round((5/9) * (curobs$TMIN - 32), 1)
    
    ## wind speeds
    curobs$WSPD <- round(curobs$WSPD * 0.44704, 1)
    curobs$WMAX <- round(curobs$WMAX * 0.44704, 1)
    
    ## rain
    curobs$RAIN <- round(curobs$RAIN * 25.4, 1)
  }
  return(curobs)
}