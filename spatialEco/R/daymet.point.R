#' @title DAYMET point values
#' @description Downloads DAYMET climate variables for specified point and timeperiod 
#' 
#' @param lat          latitude of point (decimal degrees WGS84)
#' @param long         longitude pf point (decimal degrees WGS84)
#' @param start.year   First year of data
#' @param end.year     Last year of data
#' @param site         Unique identification value that is appended to data
#' @param clean.up     (TRUE/FALSE) Delete downloaded csv files on disk
#' @param echo         (TRUE/FALSE) Echo progress
#'
#' @return A data.frame with climate results 
#'
#' @note Function uses the Single Pixel Extraction tool and returns year, yday, dayl(s), prcp (mm/day), srad (W/m^2), swe (kg/m^2), tmax (deg c), tmin (deg c), vp (Pa)
#' @note Metadata for DAYMET single pixel extraction: \url{ https://daymet.ornl.gov/files/UserGuides/current/readme_singlepointextraction.pdf }
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'         
#' @examples
#' \dontrun{
#' ( d <- daymet.point(lat = 36.0133, long = -84.2625, start.year = 2013, end.year=2014, 
#'                     site = "1", clean.up = TRUE, echo = FALSE) )
#' }
#'
#' @export
daymet.point <- function (lat, long, start.year, end.year, site, 
                          clean.up = TRUE, echo = FALSE) {
    if(missing(lat)) stop("Please define lat") 
    if(missing(long)) stop("Please define long") 
	if(missing(start.year)) stop("Please define start year") 
    if(missing(end.year)) stop("Please define end year") 
    if (start.year < 1980) stop("Data not avaliable prior to 1980")
    year.range = paste(seq(start.year, end.year, by = 1), collapse = ",")
    download.url = sprintf("https://daymet.ornl.gov/data/send/saveData?lat=%s&lon=%s&measuredParams=tmax,tmin,dayl,prcp,srad,swe,vp&year=%s", 
                           lat, long, year.range)
	if(missing(site)) site = "tmp"
    daymet.file = paste(site, "_", start.year, "_", end.year, ".csv", sep = "") 
	  if (echo == TRUE) {
          cat(paste("Downloading DAYMET data for: ", site, " at ", 
              lat, "/", long, " latitude/longitude !\n", sep = ""))
      }
    options(warn=-1)	  
	try( utils::download.file(download.url, daymet.file) )
    options(warn=0)    
      if (file.exists(daymet.file)) {
        dat = try(utils::read.table(daymet.file, header = T, sep = ",", skip = 7), silent = TRUE)
          if (inherits(dat, "try-error")) {
            unlink(daymet.file)
            stop("There is an issue with these coordinates or data request")
          }
      }
	names(dat) <- c("year", "julian", "dayl", "prcp", "srad", "swe", "tmax", "tmin", "vp")
	  if(!missing(site)) dat <- data.frame(site = site, dat)
    if(clean.up) unlink(daymet.file)	
  return(dat)
 }
