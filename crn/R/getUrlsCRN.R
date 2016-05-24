getUrlsCRN <- function(url = CRN.DAILY.URL, year = 2011){
  require("RCurl")
  yearUrl <- paste(url,as.character(year),"/",sep = "")
  ftpDir  <- getURL(yearUrl,ftp.use.epsv = FALSE, ftplistonly = TRUE)
  names <- strsplit(ftpDir,split = "\r\n") 
  allurls <- mapply(paste,yearUrl,names,sep = "")
  return(allurls)
}