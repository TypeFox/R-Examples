## COERCE HEADERS of TAXON TABLE to STANDARD
## package: megaptera
## author: Christoph Heibl
## last update: 2014-03-10

sqlTaxonomyHeader <- function(tax){
  
  ## enforce standard column names
  ## -----------------------------
  names(tax) <- tolower(names(tax))
  names(tax) <- gsub("species", "spec", names(tax))
  names(tax) <- gsub("genus|genera", "gen", names(tax))
  names(tax) <- gsub("family", "fam", names(tax))
  names(tax) <- gsub("order", "ord", names(tax))
  names(tax) <- gsub("class", "class", names(tax))
  names(tax) <- gsub("[.]", "_", names(tax))
  tax
}
