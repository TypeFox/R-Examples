#' tide.info gets the annual tide calendar for one particular location.
#' @title Annual tide calendar for one particular location
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return Return a data.frame with tide calendar:\cr
#' Level is the tide level, Tide is the High or Low Tide information and Date.Time is
#' the date/time in POSIXlt format.
#' @param file An html file from the site \code{http://tides.mobilegeographics.com/}
#' @param location Code based on \code{http://tides.mobilegeographics.com/}
#' @param year Year to get the calendar
#' @param latitude The latitude of the tide information
#' @param longitude The longitude of the tide information
#' @param tz Timezone
#' @family Periodic patterns of indices
#' @description The script extracts tide information from \cr
#' \code{http://tides.mobilegeographics.com/} into a data.frame.\cr
#' The presence of XLM package is required for this function.\cr
#' @keywords Tide
#' @examples
#' \dontrun{
#' library("HelpersMG")
#' lat <- 5.74
#' long <- -54
#' Awala2004 <- tide.info(year=2004, longitude=long, latitude=lat, tz="America/Cayenne")
#' with(Awala2004, plot(Date.Time, Level, bty="n", las=1, type="l", 
#' xlab=paste("Year", as.POSIXlt(Date.Time[1])$year+1900), 
#' ylab="Tide level in m"))
#' }
#' @export

tide.info <- function(file=NULL, year=as.POSIXlt(Sys.time())$year+1900, 
	location=0, latitude=NA, longitude=NA, tz="") {
        
  if (!requireNamespace("XML", quietly = TRUE)) {
    stop("XML package is necessary for this function")
  }
  
  # file=NULL; year=as.POSIXlt(Sys.time())$year+1900;location=0;latitude=NA;longitude=NA;tz=""
  # longitude=4.01; latitude=6.4
  if (!is.na(latitude) & !is.na(longitude)) {
    locationTide <- NULL
    load(file=paste0(system.file("Tide", package="HelpersMG"), "/locationTide.rda"))
  location <- locationTide$location[which.min(abs(locationTide$latitude-latitude)+abs(locationTide$longitude-longitude))][1]
  }
  rm(locationTide)
  if (is.null(file)) {
  theurl <- paste0("http://tides.mobilegeographics.com/calendar/year/", location, ".html?y=", year, "&m=1&d=1")
} else {
  theurl <- file
}
tables <- XML::readHTMLTable(theurl, stringsAsFactors=FALSE)
n.rows <- lapply(tables, function(t) dim(t)[1])
n.rows <- lapply(n.rows, function(t) {ifelse(is.null(t), 1, t)})
n.rows <- unlist(n.rows)

tl <- Sys.getlocale(category = "LC_TIME")
Sys.setlocale(category = "LC_TIME", locale = gsub(".._..(.+)", "en_US\\1", tl) )
Tide.Calendar <- NULL
months <- which(n.rows>27)
for (month in 1:12) {
  table <- tables[[months[month]]]
  Date <- paste0(year,"-", ifelse(month<10, "0", ""), as.character(month), "-", gsub(" ", "0", gsub("... (..)", "\\1", table[,1])))
  # Read morning high tide
  for (col in 2:6) {
  Time <- gsub("([0-9]+:[0-9]+ [AP]M) .+", "\\1", table[,col])
  
  Date.Time <- strptime(paste(Date, Time), format="%Y-%m-%d %I:%M %p", tz=tz)
  Level <- gsub("[0-9]+:[0-9]+ [AP]M .+ / (-?[\\.0-9]+) m", "\\1", table[,col])
  Tide.Calendar <- rbind(Tide.Calendar, data.frame(Date.Time=Date.Time, Level=as.numeric(Level), 
                                                   Tide=c("High Tide", "Low Tide", "High Tide", 
                                                          "Low Tide", "High Tide")[col-1], 
                                                   stringsAsFactors=FALSE))
  }
}
Tide.Calendar <- na.omit(Tide.Calendar)
Tide.Calendar <- Tide.Calendar[order(Tide.Calendar[,1]),]
Sys.setlocale(category = "LC_TIME", locale = tl )
return(Tide.Calendar)
}
