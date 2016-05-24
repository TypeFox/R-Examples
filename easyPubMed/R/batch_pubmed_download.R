batch_pubmed_download <-
function(pubmed_query_string, dest_dir = NULL, batch_size = 400, res_cn = 1){
  #handle directories
  baseDir <- getwd()
  if (!is.null(dest_dir)) {
    setwd(as.character(dest_dir))
  }
  #start by performing a PubMed ID query
  myQuery <- NULL
  while (is.null(myQuery)) {
    myQuery <- tryCatch(get_pubmed_ids(pubmed_query_string), error = function(e) NULL)
  }
  pubsNum <- as.numeric(myQuery$Count)
  tmpPapers <- NULL
  myRetstart <- 0
  myRetmax <- batch_size
  j = 1
  expTot <- pubsNum / batch_size
  if (expTot > as.integer(expTot)) {
    expTot <- as.integer(expTot) + 1
  } else {
    expTot <- as.integer(expTot)
  }
  while (myRetstart < pubsNum){
    if (j < res_cn){
      message(paste("cycle", j, "/", expTot, "skipped...", sep = " "))
    } else {
      while (is.null(myQuery) | is.null(tmpPapers)) {
        myQuery <- tryCatch(get_pubmed_ids(pubmed_query_string), error = function(e) NULL)
        tmpPapers <- tryCatch(fetch_pubmed_data(myQuery, retstart = myRetstart, retmax = myRetmax),
                              error = function(e) NULL,
                              finally = print(paste("XML batch", j , "/", expTot, "downloaded...", sep = " ")))
        if (is.null(tmpPapers)) {
          message("XML error. Retrying...")
        }
      }
      totDigits <- nchar(as.character(expTot)) + 1
      myExt <- paste(rep(0,totDigits - nchar(as.character(j))), collapse="")
      doSaveXml <- tryCatch(saveXML(tmpPapers, paste("xml_papers_", myExt, j, ".xml",sep="")), error = function(e) "ERROR")
      myQuery <- NULL
      tmpPapers <- NULL
      if (doSaveXml == "ERROR") {
        myRetstart <- myRetstart - myRetmax
        j <- j - 1
        message("An error occurred... Downloading XML again...")
      }
    }
    myRetstart <- myRetstart + myRetmax
    j <- j + 1
  }
  setwd(baseDir)
}
