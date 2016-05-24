getEmptyCsv <- function(directory = "EnvCanada", Monthly = MONTHLY.STATION.LIST){
  
  filePat <- "^99.+Env\\.csv$"
  Information <- file.info(list.files(path = directory, full.names = TRUE,pattern = filePat))
   
  dex   <- which(Information$size == 0)
  if (length(dex) == 0) {
    print("No empty files")
    return(NULL)
  }
  stationlist <- readMonthlyCsv(Monthly)
  emptyFiles <- Information[dex,]
  dirChar <- nchar(directory) + 2
  ids <- as.numeric(substring(rownames(emptyFiles),dirChar, dirChar + 7)) 
  Get <- which(stationlist$Id %in% ids)
  return(Get)
  
}