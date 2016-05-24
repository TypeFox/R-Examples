fetch_pubmed_data <-
function(pubmed_id_list, retstart = 0, retmax = 500) {
  myIDlist <- pubmed_id_list
  # quality checks
  if ((!is.list(myIDlist)) |
      is.na(myIDlist$WebEnv) |
      is.na(myIDlist$QueryKey) |
      is.na(myIDlist$Count) |
      !is.integer(as.integer(retstart)) |
      !is.integer(as.integer(retmax))
  ) {
    stop("There is an issue with the PubMed ID list you supplied. Please, call the function again and supply the result of a <get_pubmed_ids()> call as argument. Thank you.")
  } else {
    
    # retrieve what we need from the Pubmed ID List
    myWebEnv <- myIDlist$WebEnv
    myKey <- myIDlist$QueryKey
    myCount <- as.numeric(as.character(myIDlist$Count))
    
    # retrieve what we need from the function call (function arguments parsing)
    myRetstart = as.integer(retstart)
    if (myRetstart < 0) { myRetstart = 0 }
    myRetmax <- as.integer(retmax)
    # limit retmax to 5000
    if (myRetmax > 5000) { myRetmax = 5000 }
    if (myRetmax < 1) { myRetmax = 1 }
    
    # put together the eutils URL for the query -- will retrieve a XML result
    efetch_url = paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
                       "db=pubmed&WebEnv=",myWebEnv,
                       "&query_key=",myKey,
                       "&retstart=",myRetstart,
                       "&retmax=",myRetmax,
                       "&rettype=null&retmode=xml",sep="")
    
    # connect and fetch the XML data - no files are saved. We do everything on the fly
    tmpConnect <- url(efetch_url, open = "rb")
    abstrXML <-readLines(tmpConnect, warn = FALSE, encoding = "xml")
    abstrXML <- xmlParse(abstrXML)
    close.connection(tmpConnect)
    
    # retrun the XML data we retrieved from Pubmed. It contains all the data
    return(abstrXML)
  }
}
