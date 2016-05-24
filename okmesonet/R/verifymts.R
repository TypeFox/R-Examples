verifymts<- function(station, datemts) {
  ## Verify access to MTS file
  ##
  ## Arguments:
  ##  station: four letter character ID for Mesonet station, lowercase
  ##  datemts: date of MTS file to retrieve, POSIXct format
  ## Returns: logical indicating success
  mtspath <- paste("http://www.mesonet.org/index.php/dataMdfMts/",
                   "dataController/getFile/", 
                   format.POSIXct(datemts, format="%Y%m%d"), station, 
                   "/mts/TEXT/", sep = "")
  internetoption <- getOption("internet.info")
  options(show.error.messages=FALSE, internet.info=3)
  checkconn <- try(readLines(mtspath, n=1))
  options(show.error.messages=TRUE, internet.info=internetoption)
  if(inherits(checkconn,"try-error")) return(FALSE) else return(TRUE)
}