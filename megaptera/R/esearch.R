## package: megaptera
## called by: 
## author: Christoph Heibl (at gmx.net)
## last update: 2014-08-20

esearch <- function(tax, db = "taxonomy", usehistory = FALSE, retmax = 10000){
  
  
  id <- seq(from = 1, to = length(tax), by = 100)
  id <- data.frame(from = id, to = c(id[-1] - 1, length(tax)))
  webEnv <- NULL; queryKey <- NULL
  
  for (i in 1:nrow(id) ){
    
    cat("\niteration: ", i)
    term <- unlist(tax)[id[i, 1]:id[i, 2]]
    term <- paste(term, collapse = "+OR+")
    
    ## assemble URL
    xml <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/", 
                 "eutils/esearch.fcgi?",
                 "tool=megaptera",
                 "&email=heibsta@gmx.net&",
                 "&usehistory=y&", sep = "")
    if ( !is.null(webEnv) ) xml <- paste(xml, "&WebKey=", webEnv, sep = "")
    xml <- paste(xml, "&db=", db,  "&term=", term, 
                 "&retmax=", retmax, sep = "")
    
    ## get and parse result
    if ( is.null(webEnv) ){
      xml <- xmlTreeParse(xml, useInternalNodes = TRUE)
      webEnv <- xpathSApply(xml, fun = xmlToList,
                            path = "//eSearchResult/WebEnv")
      queryKey <- xpathSApply(xml, fun = xmlToList,
                              path = "//eSearchResult/QueryKey")
    } else {
      browseURL(xml, browser = "false")
    }  
  }
  
#   s <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
#              "esummary.fcgi?tool=megaptera&email=heibsta@gmx.net&",
#              "db=nucleotide&query_key=", queryKey, 
#              "&WebEnv=", webEnv, sep = "")
#   s <- xmlTreeParse(s, useInternalNodes = TRUE)
#   ss <- xpathSApply(s, fun = xmlToList,
#                     path = "//DocSum/Id")
  
  x <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
             "efetch.fcgi?tool=megaptera&email=heibsta@gmx.net&",
             "db=nucleotide&query_key=", queryKey, 
             "&WebEnv=", webEnv,
             "&rettype=gb&retmode=xml", sep = "")
  #   xml <- scan(x, what = "c", quiet = TRUE, sep = "\n")
  #   write(xml, "test2.xml"); system("open -t test2.xml")
  #   x  <- "test.xml"
  obj <- xmlTreeParse(x, getDTD = FALSE, useInternalNodes = TRUE)
  obj
}



