## package: megaptera
## called by: stepA
## author: Christoph Heibl (at gmx.net)
## last update: 2014-10-25

term <- function(organism, locus) {
  
  if ( !is.character(organism) )
    stop("organism must be of mode 'character'")
  if ( !inherits(locus, "locus") )
    stop("locus must be of class 'locus'")
  
  organism <- paste(organism, "[orgn]", sep = "")
  ## eventually move to higher level:
  aliases <- paste("\"", locus@aliases, "\"", sep = "")
  
  ## create URL using 'sgene' object
  ## -------------------------------
  URL <- vector(length = length(organism))
  for ( i in seq_along(organism) ){
    url <- paste(organism[i], "AND",
                 NCBI.wrap(aliases, field = locus@search.fields))
    if ( !"not" %in% locus@not ){
      not <- paste("\"", locus@not, "\"", sep = "")
      url <- paste(url, "NOT", 
                   NCBI.wrap(not, field = locus@search.fields))
    }
    if ( length(url) > 1 ){
      url <- paste("(", url, ")", sep = "")
      url <- paste(url, collapse = " OR ")
    }
    URL[i] <- url
  } # end of FOR-loop over i
  if ( length(URL) > 1 ){
    URL <- paste("(", URL, ")", sep = "")
    URL <- paste(URL, collapse = " OR ")
  }
  URL
}
