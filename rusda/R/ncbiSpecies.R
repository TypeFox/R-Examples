##                       ncbiSpecies                 ##
##      This code is part of the rusda package       ##
##     F.-S. Krah 2016 (last update: 2016-03-18)     ##

ncbiSpecies <- function(sciname, clean, sub){
  uid <- get_uid(sciname, check = TRUE)[[1]]
  url <- paste("http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=",uid, sep="")
  parse <- XML::htmlTreeParse(url, useInternal = TRUE)
  parse <- xpathApply(parse, "//strong", xmlValue)
  sp <- unlist(parse)
  sp <- sp[grep(sciname, sp)]
  if(clean == TRUE)
  {
    nowant <- c("f\\.|sp\\.|var\\.|subsp\\.|ssp\\.|x|unclassified")
    if(length(grep(nowant, sp))>0)
    {sp <- sp[-grep(nowant, sp)]}
  }
  if(!sub == TRUE)
  {
    sp <- sp[vapply(strsplit(sp, "\\W+"), length, integer(1))==2]
  }
  return(sp)
}
