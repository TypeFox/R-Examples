getMissingScrape <- function(monthlyList = MONTHLY.STATION.LIST, directory = "EnvCanada"){
  if (!file.exists(directory)) stop(" Directory does not exist")
  monthly <- read.csv(monthlyList)
  allId   <- monthly$Id  
  filePat <- "^99.+Env\\.csv$"
  csvFiles <- list.files(path = directory, full.names = TRUE, pattern = filePat)
  stations <- length(csvFiles)
  dirChar <- nchar(directory) + 2
  ids <- as.numeric(substring(csvFiles,dirChar, dirChar + 7))
  missing <- setdiff(allId, ids)
  dex <- which(allId %in% missing)
  if (length(dex) == 0) dex <- NULL
  return(dex)
  
  
}