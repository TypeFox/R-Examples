createInventory   <- function(directory = "EnvCanada"){
  
  if (!file.exists(directory)) stop(" Directory does not exist")
  filePat <- "^99.+Env\\.csv$"
  csvFiles <- list.files(path = directory,full.names = TRUE, pattern = filePat)
  stations <- length(csvFiles)
  dirChar <- nchar(directory)+2
  ids <- as.numeric(substring(csvFiles,dirChar,dirChar+7))
  Inventory <- data.frame(Name      = rep(NA,stations),
                          Province  = rep(NA,stations),
                          Lat       = vector(mode = "numeric",   length = stations),
                          Lon       = vector(mode = "numeric",   length = stations),
                          Altitude  = vector(mode = "numeric",   length = stations),
                          Id        = rep(NA,stations),
                          WMO       = vector(mode = "numeric",   length = stations),
                          TCid      = rep(NA,stations),
                          downloadId        = ids) 
 
     
    for( files in 1 :length(csvFiles)){
  
        print(csvFiles[files]) 
        S      <- read.csv(csvFiles[files], nrows = 7, stringsAsFactors = FALSE)
        inv    <- as.data.frame(t(S),stringsAsFactors=FALSE)
   
        Inventory$Name[files]      <- rownames(inv)[2]
        Inventory$Province[files]  <- inv$V1[2]
        Inventory$Lat[files]       <- inv$V2[2]
        Inventory$Lon[files]       <- inv$V3[2]
        Inventory$Altitude[files]  <- inv$V4[2]
        Inventory$Id[files]        <- inv$V5[2]
        Inventory$WMO[files]       <- inv$V6[2]
        Inventory$TCid[files]      <- inv$V7[2]
         
       
   }
  return(Inventory)
} 
  
 