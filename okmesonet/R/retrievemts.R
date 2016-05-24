retrievemts <- function(station, datemts, variables, missingNA) {
  ## Retrieve MTS file from Mesonet website
  ## Arguments:
  ##  station: four letter character ID for Mesonet station, lowercase
  ##  datemts: date of MTS file to retrieve, POSIXct format
  ##  variables: variables to retrieve
  ## Returns: dataframe containing MTS file with timestamp
  
  ## read MTS from Mesonet website
  mtspath <- paste("http://www.mesonet.org/index.php/dataMdfMts/",
                   "dataController/getFile/", 
                   format.POSIXct(datemts, format="%Y%m%d"), station, 
                   "/mts/TEXT/", sep = "")
  mts <- read.csv(mtspath, skip = 2, header = T, as.is = T, sep = "", 
                  nrows = 288)
  
  ## IMPORTANT: convert 'TIME' field to timestamp
  ## TIME represents the number of minutes from base time specific 
  ## in MTS file
  ## see http://www.mesonet.org/wiki/Public:MDF_Format
  ## this appears to be always 00:00:00 UTC
  mts$TIME <- as.POSIXct(mts$TIME*60, origin=datemts)
  
  ## reorder data frame with STID, STNM, and timestamp first  
  ## variables to put at end
  endvar <- colnames(mts)[names(mts)!="STID" & names(mts)!="STNM" 
                           & names(mts)!="TIME"]
  
  ## replace data frame with desired column order by subsetting
  mts <- mts[, c("STID", "STNM", "TIME", endvar)]
  
  ## replace values <=-994 with NA
  if (missingNA == T) mts <- replace(mts, mts<=c(-994), NA)
  
  ## return data frame with desired vairables
  if(all(variables %in% "ALL")) {
    return(mts)
  } else {
    return(mts[, c("STID", "STNM", "TIME", variables)])
  }
}