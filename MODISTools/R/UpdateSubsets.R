UpdateSubsets <- 
function(LoadDat, StartDate = FALSE, Dir = ".")
{
	if(StartDate) details <- LoadDat[!duplicated(data.frame(LoadDat$lat, LoadDat$long, LoadDat$end.date, LoadDat$start.date)), ]
  if(!StartDate) details <- LoadDat[!duplicated(data.frame(LoadDat$lat, LoadDat$long, LoadDat$end.date)), ]
	cat("Found", nrow(details), "unique time-series in original file\n") 
	
	# Year or posixt date format?
	Year <- FALSE
	POSIXt <- FALSE
	posix.compatible <- try(as.POSIXlt(details$end.date), silent = TRUE)
	if(any(class(details$end.date) == "POSIXt") | all(class(posix.compatible) != "try-error")) POSIXt <- TRUE
	if(all(is.numeric(details$end.date) & nchar(details$end.date) == 4) & 
	     any(class(posix.compatible) == "try-error")) Year <- TRUE
	if(!Year & !POSIXt) stop("Date information in LoadDat is not recognised as years or as POSIXt format.")
	if(Year & POSIXt) stop("Date information in LoadDat is recognised as both year and POSIXt formats.")
  
	if(Year) endyear <- details$end.date
	if(POSIXt) endyear <- as.numeric(format(details$end.date, "%Y"))
	if(StartDate){
		if(Year) startyear <- details$start.date
		if(POSIXt) startyear <- as.numeric(format(details$start.date, "%Y"))
	}
	
	ID <- ifelse(any(names(details) == "ID"), TRUE, FALSE)
	fmt <- '%.5f'
	if(StartDate){
    if(ID){
      ## Check that all author-given IDs will be unique for each unique time-series, and check that they won't cause issues with product information
      n.unique <- length(unique(details$ID)) == nrow(details)
      if(!n.unique){
        cat('Number of IDs is not unique.\n')
        details$ID <- paste("Lat", sprintf(fmt, details$lat), "Lon", sprintf(fmt, details$long), "Start", startyear, "End", endyear, sep = "")
      }
    } else {
      details$ID <- paste("Lat", sprintf(fmt, details$lat), "Lon", sprintf(fmt, details$long), "Start", startyear, "End", endyear, sep = "")
    }
	}
	if(!StartDate){
    if(ID){
      ## Check that all author-given IDs will be unique for each unique time-series, and check that they won't cause issues with product information
      n.unique <- length(unique(details$ID)) == nrow(details)
      if(!n.unique){
        details$ID <- paste("Lat", sprintf(fmt, details$lat), "Lon", sprintf(fmt, details$long), "End", endyear, sep = "")
      }
    } else {
      details$ID <- paste("Lat", sprintf(fmt, details$lat), "Lon", sprintf(fmt, details$long), "End", endyear, sep = "")
		} 	
	}
	
	filelist <- list.files(path = Dir, pattern = ".asc")
	cat("Found", length(filelist), "subsets previously downloaded\n")
	
	char.position <- regexpr("___", filelist)
	downloaded.ids <- substr(filelist, 1, char.position -1)
	
	revised.subsets <- details[which(details$ID %in% downloaded.ids == FALSE), ]
	return(revised.subsets)
}