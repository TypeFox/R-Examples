#' Standardize taxonomic names
#' 
#' This function standardizes taxa names. It is used mainly internally, but might be
#' helpful to the end user in some situations.
#' 
#' @param taxon a character vector containing a single name
#' @return a character vector
#' @export
#' @examples
#' \dontrun{
#' standardize.names("Miconia sp 01")
#' standardize.names("Miconia Sp 2")
#' standardize.names("Sp18")
#' }

standardize.names <- function(taxon) {
  #taxon <- iconv(taxon, to = "ASCII//TRANSLIT")
  taxon <- gsub("[I|i]ndet[a-z]*", "Indet\\.", taxon)
  isolated <- grepl("^[s|S][P|p]\\.+", taxon)
  if (isolated) {
    taxon <- paste("Indet.", taxon)
  }
  taxon.split <- strsplit(taxon, "\\s+")[[1]]
  taxon.number <- regmatches(taxon, regexpr("[0-9]+", taxon))
 
  if (length(taxon.number) == 0L) {
    taxon.ident <- regmatches(taxon.split[2], regexpr("[s|S][P|p]\\.*$", taxon.split[2]))
    if (length(taxon.ident) == 0L) {
      return(taxon)
    } else {
      taxon.std <- paste(taxon.split[1], "sp.")
      return(taxon.std)
    }
  } else {
    taxon.number <- as.numeric(taxon.number)
    taxon.std <- paste(taxon.split[1], " sp.", taxon.number, sep = "")
    return(taxon.std)
  }
}