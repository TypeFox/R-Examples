get_pubmed_ids <-
function(pubmed_query_string) {
  # coerce PubMed Query to character vector(String)
  myQuery <- as.character(pubmed_query_string)
  myQuery <- gsub(" ", "+", myQuery, fixed = TRUE)
  # Pubmed API - query URL
  myPubmedURL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
                       "db=pubmed&term=",
                       myQuery,
                       "&usehistory=y",sep='')
  # connect to Pubmed and fetch the PubMed IDs corresponding to the Query
  IDconnect <- url(myPubmedURL, open = "rb")
  idXML<-readLines(IDconnect, warn = FALSE, encoding = "xml")
  idXML<- xmlTreeParse(idXML)
  close.connection(IDconnect)
  myIDlist <- xmlToList(idXML)
  return(myIDlist)
}
