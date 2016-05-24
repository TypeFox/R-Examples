setClass("locus", 
         representation = list(
           sql = "character", 
           aliases = "character", 
           not = "character",
           search.fields = "character")
)

"locus" <- function(..., not,
                    search.fields = c("gene", "title"),
                    check = TRUE){
  if ( missing(...) ){
    new("locus", 
        aliases = "undefined", 
        not = "undefined", 
        sql = "undefined",
        search.fields = search.fields
    )
  } else {
    
    aliases <- c(...)
    
    ## GenBank uses uppercase and lowercase spelling
    ## in the same places ...
    aliases <- unique(c(aliases, 
                        toupper(aliases),
                        tolower(aliases), 
                        paste(toupper(substring(aliases, 1, 1)), 
                              substring(aliases, 2), sep = "")))
    
    if ( missing(not) ) not <- "not"
    
    ## check if locus exists at GenBank
    ## --------------------------------
    if ( check ){
      cat("\nChecking if locus exists on GenBank ..")
      a <- paste("\"", aliases, "\"", sep = "")
      url <- paste(aliases, "[all]", collapse = " OR ")
      if ( !"not" %in% not ) {
        n <- paste("\"", not, "\"", sep = "")
        url <- paste(url, "NOT", NCBI.wrap(not, field = "all"))
      }
      url <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
                   "esearch.fcgi?db=nucleotide&term=", url, 
                   "&rettype=gb&retmode=xml", sep = "")
      x <- xmlTreeParse(url, getDTD = FALSE, 
                        useInternalNodes = TRUE)
      x <- unique(xpathSApply(x, "//Count", xmlValue))
      cat("\n.. found", x[1], "records")
    }
    new("locus", 
        aliases = aliases, 
        not = not, 
        sql = sql.conform(aliases[1]),
        search.fields = search.fields
    )
  }
}

setMethod("show",
          signature(object = "locus"),
          function (object) 
          {
            if ( object@sql == "undefined" ){
              cat("\nLocus definition: empty")
            } else {
              cat("\nLocus definition: ", object@aliases[1],
                  "\nThe strings", 
                  paste("\n -", object@aliases[-1]),
                  "\nwill be searched for in the search fields", 
                  paste("\n -", object@search.fields),
                  "\nSequences will be stored in tables '", 
                  paste("acc", object@sql, sep = "_"),
                  "' and '",
                  paste("spec", object@sql, sep = "_"),
                  "'", sep = "")
            }
          }
)
