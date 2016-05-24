# evaluate whether the URL is available
is.URLconnectable <- function(theURL) {
  aConnection <- suppressWarnings(try(url(theURL, open = "rb"), silent = TRUE))
  if(inherits(aConnection, "try-error")) return(FALSE)
  close(aConnection)
  return(TRUE)
}

getHTMLfromURL <- function(theURL) {
  theHTMLvector <- tryCatch(
    readLines(theURL, n = -1, ok = TRUE, warn = FALSE), 
    error = function(cond) return(paste0(" failed url, with error: ", cond[1])),
    warning = function(cond) return(paste0(" maybe failed url, with warning: ", cond[1]))
  )    
  return(theHTMLvector)
}

getPDFfromURL <- function(theURL, 
                          theDirectory, 
                          theFileName) {
  aConnection <- suppressWarnings(try(download.file(theURL, 
                                                    destfile = file.path(theDirectory, 
                                                               paste0(theFileName, ".pdf")), 
                                                    quiet = TRUE, 
                                                    method = "auto",  
                                                    mode = "wb", 
                                                    cacheOK = FALSE), silent = TRUE))
  if(inherits(aConnection, "try-error")) 
    return(paste0(" failed download, with error: ", aConnection[1]))
  
  return(aConnection)
}

extractPDFsFromHTML <- function(theHTMLvector, 
                                theDirectory, 
                                theFileName, 
                                validatePDF = TRUE) {
  
  candidateLinks <- extractPDFLinksFromHTML(theHTMLvector, validatePDF)
  if(candidateLinks[1] == FALSE) return (" failed, no valid url links detected")
  
  if(length(candidateLinks) != 1) theFileName <- paste0(theFileName, 
                                                        ".", 
                                                        1:length(candidateLinks))
  
  downloadOutcomes <- mapply(
    getPDFfromURL, 
    candidateLinks, 
    theDirectory, 
    theFileName
  )
  
  #return TRUE if at least one file was successfully downloaded
  return(any(downloadOutcomes == 0))
}

extractPDFLinksFromHTML <- function(theHTMLvector, 
                                    validatePDF = TRUE) {
  
  # extract all unique PDF-links using a collection of search criteria defined in wildcardFunctionList
  candidateLinks <- unique(unlist(lapply(wildcardFunctionList, 
                                         function(x) eval(parse(text = paste0(x, "(theHTMLvector)"))))))

  # from the candidateLinks, select only those with connectable URLs
  candidateLinks <- candidateLinks[unlist(lapply(candidateLinks, is.URLconnectable))]
  
  # from the candidateLinks, select only those identified as PDFs
  if(validatePDF) candidateLinks <- candidateLinks[unlist(lapply(candidateLinks, isPDF))]

  if(length(candidateLinks) == 0) return (FALSE)
  return(candidateLinks)
}

# collection of HTML PDF-link extractor functions
wildcardFunctionList <- c(
  "wcf_generic",  # generic links
  "wcf_jstore",   # jstore links (www.jstor.org)
  "wcf_bioone",   # bioone links (www.bioone.org)
  "wcf_wiley",    # wiley online links (onlinelibrary.wiley.com) NOTE: no longer works as of 11/5/14
  "wcf_nrc",      # NRC publications links (http://www.nrcresearchpress.com)
  "wcf_jstage",   # J-stage publications links (www.jstage.jst.go.jp)
  "wcf_tfo",      # Taylor & Francis online links (www.tandonline.com)
  "wcf_nature",   # nature publications links (www.nature.com) 
  "wcf_maney",    # Maney online (www.maneyonline.com)
  "wcf_acs"       # acs publications links (pubs.acs.org)
)

# extract all generic PDF-links from HTML vector document
# last validated: 11/3/2014, M.J. Lajeunesse
wcf_generic <- function(theHTMLdata) {
  someLinks <- pullLinks(theHTMLdata, "http.*pdf")
  candidateLinks <- vector()
  for(i in 1:length(someLinks)) {
    cleanedURL <- wcf_generic_URL_cleaner(someLinks[i])
    if(!identical(cleanedURL, character(0))) candidateLinks <- c(candidateLinks, cleanedURL)
  }
  if(length(candidateLinks) == 0) return(NULL)
  else return(candidateLinks)
}

# HELPER FUNCTION for wcf_generic(): extra scrubbing of generic html links 
wcf_generic_URL_cleaner <- function(aURLstring) {
  subStrings <- unlist(strsplit(aURLstring, " "))
  findURLS <- str_extract(subStrings, "http.*\\.pdf")
  return(unique(findURLS[!is.na(findURLS)]))
}

# extract all PDF-links from a HTML vector document from jstore
# last validated: 11/3/2014, M.J. Lajeunesse
wcf_jstore <- function(theHTMLdata) {
  candidateLinks <- pullLinks(theHTMLdata, "pdfplus.*pdf")
  if(length(candidateLinks) == 0) return(NULL)
  candidateLinks <- paste0("http://www.jstor.org/stable/", 
                           candidateLinks, "?acceptTC=true")
  return(candidateLinks)
}

# extract all PDF-links from a HTML vector document from bioone
# last validated: 11/3/2014, M.J. Lajeunesse
wcf_bioone <- function(theHTMLdata) {
  candidateLinks <- pullLinks(theHTMLdata, "/doi/pdf/.*\">PDF")
  if(length(candidateLinks) == 0) return(NULL)
  candidateLinks <- paste0("http://www.bioone.org",  
                           gsub("\">PDF", "", candidateLinks))
  return(candidateLinks)
}

# extract all PDF-links from a HTML vector document from wiley online
# last validated: no longer works as of 11/3/2014, M.J. Lajeunesse
# NOTE: can download pdf through browseURL(candidateLinks), but too much work to
# locate, move, and rename downloaded pdf
#wcf_wiley <- function(theHTMLdata) {
#  candidateLinks <- gsub("abstract", "pdf", unique(pullLinks(theHTML, 
# "http://onlinelibrary.wiley.com/doi/.*/abstract")))
#  return(candidateLinks)
#}

# extract all PDF-links from a HTML vector document from ACS publications
# last validated: 11/3/2014, M.J. Lajeunesse
wcf_acs <- function(theHTMLdata) {
  candidateLinks <- pullLinks(theHTMLdata, "/doi/pdf/.*\">PDF<")
  if(length(candidateLinks) == 0) return(NULL)
  candidateLinks <- paste0("http://pubs.acs.org", 
                           gsub("\">PDF<", "", candidateLinks))
  return(candidateLinks)
}

# extract all PDF-links from a HTML vector document from NRC publications
# last validated: 11/3/2014, M.J. Lajeunesse
wcf_nrc <- function(theHTMLdata) {
  candidateLinks <- pullLinks(theHTMLdata, "/doi/pdf/.*\">")
  if(length(candidateLinks) == 0) return(NULL)
  candidateLinks <- paste0("http://www.nrcresearchpress.com", 
                           gsub("\">", "", candidateLinks))
  return(candidateLinks)
}

# extract all PDF-links from a HTML vector document from J-stage publications
# last validated: 11/3/2014, M.J. Lajeunesse
wcf_jstage <- function(theHTMLdata) {
  candidateLinks <- pullLinks(theHTMLdata, "/article.*/_pdf")
  if(length(candidateLinks) == 0) return(NULL)
  candidateLinks <- paste0("https://www.jstage.jst.go.jp", candidateLinks)
  return(candidateLinks)
}

# extract all PDF-links from a HTML vector document from Taylor & Francis online publications
# last validated: 11/3/2014, M.J. Lajeunesse
wcf_tfo <- function(theHTMLdata) {
  candidateLinks <- pullLinks(theHTMLdata, "/doi/pdf/.*\">D")
  if(length(candidateLinks) == 0) return(NULL)
  candidateLinks <- paste0("http://www.tandfonline.com", 
                           gsub("\">D", "", candidateLinks))
  return(candidateLinks)
}

# extract all PDF-links from a HTML vector document from Nature publications
# last validated: 11/3/2014, M.J. Lajeunesse
wcf_nature <- function(theHTMLdata) {
  candidateLinks <- unique(unlist(strsplit(pullLinks(theHTMLdata, 
                                                     "=.*\\.pdf"), "\"")))
  if(length(candidateLinks) == 0) return(NULL)
  candidateLinks <- unique(pullLinks(candidateLinks, ".*\\.pdf"))
  candidateLinks <- paste0("http://www.nature.com", candidateLinks)
  return(candidateLinks)
}

# extract all PDF-links from a HTML vector document from Maney publications
# last validated: 11/3/2014, M.J. Lajeunesse
wcf_maney <- function(theHTMLdata) {
  candidateLinks <- pullLinks(theHTMLdata, "/.*pdfplus.*\"")
  if(length(candidateLinks) == 0) return(NULL)
  candidateLinks <- paste0("http://www.maneyonline.com", 
                           gsub("\"", "", candidateLinks))
  return(candidateLinks)
}

#wcf_wiley <- function(theHTMLdata_redirect) {
#  theHTMLdata <- getHTMLfromURL(wcf_wiley_helper(unlist(strsplit(theHTMLdata_redirect, ">"))))
#  candidateLinks <- pullLinks(unlist(strsplit(theHTMLdata,"\"")), "http://onlinelibrary.wiley.com/store/.*pdf.*")
#  #if(length(candidateLinks) == 0) return(NULL)
#  return(candidateLinks) 
#}

wcf_wiley <- function(theHTMLdata) {
  redirectLinks <- gsub("abstract", 
                        "pdf", 
                        unique(pullLinks(theHTMLdata, 
                                         "http://onlinelibrary.wiley.com/doi/.*/abstract")))
  if(length(redirectLinks) == 0) return(NULL)
  theHTMLdata <- getHTMLfromURL(redirectLinks)
  candidateLinks <- pullLinks(unlist(strsplit(theHTMLdata,"\"")), 
                              "http://onlinelibrary.wiley.com/store/.*pdf.*")
  if(length(candidateLinks) == 0) return(NULL)
  return(candidateLinks)
}



# common HTML extractor
pullLinks <- function(theHTMLdata, wildcard) {
  someLinks <- str_extract(theHTMLdata, wildcard)
  return(unique(someLinks[!is.na(someLinks)]))
}

# common HTML extractor
pullLinks2 <- function(theHTMLdata, wildcard) {
  someLinks <- grep(wildcard, theHTMLdata, value = TRUE)
  return(someLinks[!is.na(someLinks)])
}
