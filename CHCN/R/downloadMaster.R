downloadMaster <- function(url = STATION.URL, localFile = MASTER.STATION.LIST){
  require(RCurl)
   
   X<-getURLContent(url, ssl.verifypeer = FALSE,verbose = TRUE)
   Xcsv <- read.csv(textConnection(X),stringsAsFactors = FALSE)
   write.csv(Xcsv,localFile)
  
}