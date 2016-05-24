tisFromCsv <- function(csvFile,
                       dateCol = "date",
                       dateFormat = "%Y%m%d",
                       tz = "", 
                       tif = NULL,
                       defaultTif = "business",
                       save = F,
                       envir = parent.frame(),
                       naNumber = NULL,
                       chopNAs = TRUE,
                       tolerance = sqrt(.Machine$double.eps),
                       ...){
  ## csvFile is the path to a csv file with column names across the top.
  ## Everything but the first row (the column names) should be numeric, and there
  ## must as many column names (enclosed in quotes) as there are columns.
  ## The column named by dateCol should contain dates in the dateFormat format.
  ##
  ## The function reads in all the data and returns a list of named tis
  ## (Time Indexed Series) objects.  The names of the list are the column
  ## names from the csvFile.
  ## If chopNAs is TRUE, the naWindow() function is applied to each series to chop off
  ## leading and trailing NA observations, so the series in the list may not all start and
  ## end on the same dates. 
  retList <- list()
  zdf <- read.csv(csvFile, as.is = T, ...)
  zNames <- names(zdf)
  dateColIndex <- match(tolower(dateCol), tolower(zNames), nomatch = NA)
  if(is.na(dateColIndex))
    stop(paste(csvFile, "does not have a column named", dateCol, "\n"))
  if(all(toupper(zNames) == zNames)) names(zdf) <- tolower(zNames)
  
  zDateStrings <- as.character(zdf[[dateColIndex]])
  zdf <- zdf[-dateColIndex]
  z <- as.matrix(zdf[,unlist(lapply(zdf, is.numeric)), drop = F])

  if(!is.null(naNumber)){
    naSpots <- (1:length(z))[abs(z - naNumber) <= tolerance]
    naSpots <- naSpots[!is.na(naSpots)]
    z[naSpots] <- NA
  }
  cn <- colnames(z) 
  if(NCOL(z) == 0) stop("No non-NA values in file")

  if(tolower(dateFormat) == "excel")
    dateTimes <- POSIXct(jul(as.ssDate(zdf[[dateColIndex]])))
  else {
    if(length(grep("%d", dateFormat)) == 0){
      dateFormat <- paste(dateFormat, "%d")
      zDateStrings <- paste(zDateStrings, "1")
    }
    dateTimes <- as.POSIXct(strptime(zDateStrings, format = dateFormat, tz = tz))
  }
  naTimes <- is.na(dateTimes)
  if(any(naTimes)){
    warning("Some rows ignored due to NA dateTimes resulting from ",
            "attempts to convert the strings\n  ",
            paste(zDateStrings[naTimes], collapse = ", "),
            "\nto dates")
    dateTimes <- dateTimes[!naTimes]
    z <- z[!naTimes,,drop = FALSE]
  }
  if(!is.null(tif)) dtTi <- ti(dateTimes, tif = tif)
  else              dtTi <- inferTi(dateTimes)
  zStart <- dtTi[1]
  zEnd   <- tail(dtTi, 1)
  zSeries <- tis(matrix(NA, zEnd - zStart + 1, ncol(z)), start = zStart)
  zSeries[dtTi,] <- z
  class(zSeries) <- "tis"
  colnames(zSeries) <- colnames(z)
  
  if(chopNAs)
    retList <- lapply(columns(zSeries), naWindow)
  else
    retList <- columns(zSeries)
  
  if(save) assignList(retList, envir = envir)
  gc()
  
  if(save) invisible(retList)
  else     retList
}
