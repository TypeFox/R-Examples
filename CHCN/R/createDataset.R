createDataset   <- function(Ids = NULL,filename = MONTHLY.STATION.LIST,directory = "EnvCanada"){
  
  if (!file.exists(directory)) stop(" Directory does not exist")
  x <- read.csv(filename)
  
  All <- x$Id
  if ( is.null(Ids)) {
    getIds <- All
  } else {
    getIds <- Ids
  }
   
  getIds <- intersect(All,getIds)
   
  if (length(getIds) == 0) stop(" get must be a subset of the total number of stations")
  
  filePat  <- "^99.+Env\\.csv$"
  csvFiles <- list.files(path = directory,full.names = TRUE, pattern = filePat)
  stations <- length(csvFiles)
  dirChar  <- nchar(directory)+2
  ids      <- as.numeric(substring(csvFiles,dirChar,dirChar+7))
  dex      <- which(ids %in% getIds)
   
  
  BB <-data.frame()
     
   for( files in dex){
      
      Data <- read.csv(csvFiles[files], skip = 17, stringsAsFactors = FALSE, na.strings = "")
      Data <- cbind(Id = All[files],Data,Filename = csvFiles[files])
      BB <- rbind(BB,Data)
      print(csvFiles[files])
         
       
   }
  return(BB) 
} 
    
  
  
     