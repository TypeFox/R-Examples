collateHourly <- function(directory = HOURLY_DIR){
  dateTime <- as.character(Sys.Date())
  fname    <-  paste("CRN_Hourly_", dateTime, ".dat", sep = "")
  invname  <-  sub("dat","inv",fname)
   
  filenames <- list.files(path = directory, full.names = TRUE, pattern = ".txt")
  files     <- basename(filenames)
  directory <- data.frame(years = substr(files,10,13),stations = substr(files,15,100),stringsAsFactors=FALSE)
  StationNames <- unique(directory$stations)
  Metadata     <- data.frame(WBANNO = rep(NA,length(StationNames)),Lat = rep(NA,length(StationNames)), 
                              Lon = rep(NA,length(StationNames)),
                              Name = rep(NA,length(StationNames)))
                              
                              
  
  for ( thisStation in 1:length(StationNames)){
     
    dex <- which(directory$stations == StationNames[thisStation])
    fileList <- filenames[dex]
    for (thisFile in 1: length(fileList)){
      X <- read.table(fileList[thisFile], na.strings = -9999.0 )
      X <- apply(X,MARGIN = 2, FUN = function(x) ifelse(x == -9999,NA,x)) 
       
      write.table(X,file = fname,col.names = FALSE,row.names = FALSE,append = TRUE)
      Metadata$WBANNO[thisStation] <- X[1,1]
      Metadata$Lat[thisStation] <- X[1,8]
      Metadata$Lon[thisStation] <- X[1,7]
      Metadata$Name[thisStation] <- sub(".txt","",StationNames[thisStation])
    }
  }
  write.table(Metadata,file = invname,col.names = TRUE,row.names = FALSE, quote = FALSE, sep = ",")
}