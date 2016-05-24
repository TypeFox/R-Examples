#' Remove the author(s) from a taxon name.
#' 
#' This attempts to remove the authorities of a taxonomic name.
#' 
#' @param taxon a character vector containing a single taxon
#' @export
#' @return a character vector
#' @examples
#' \dontrun{
#' remove.authors("Coffea arabica L.")
#' remove.authors("Chrysophyllum argenteum subsp. nitidum (G.F.W.Meyer) T.D.Penn.")
#' }
remove.authors <- function(taxon) {
  #taxon <- iconv(taxon, to = "ASCII//TRANSLIT"
  taxon <- trim(taxon)
  ident <- regmatches(taxon, regexpr("\\s+sp\\.+\\w*", taxon))
  subsp <- regmatches(taxon, regexpr("subsp\\.|var\\.", taxon))
  if (taxon == "" | is.na(taxon)) return(NA)
  taxon.split <- strsplit(taxon, " ")
  taxon.split.unlist <- unlist(taxon.split)
  if (length(taxon.split.unlist) == 1) return(taxon)
  if (length(subsp) != 0L) {
    taxon.split.subsp <- unlist(strsplit(taxon, "\\ssubsp\\.\\s|\\svar\\.\\s"))
    taxon.split.subsp.1 <- unlist(strsplit(taxon.split.subsp[1], " "))
    taxon.split.subsp.2 <- unlist(strsplit(taxon.split.subsp[2], " "))
    res <- paste(paste(taxon.split.subsp.1[c(1, 2)], collapse = " "), subsp, taxon.split.subsp.2[1], collapse = " ")
    return(res)
  }
  if (grepl("^[a-z].*[a-z]+$", taxon.split.unlist[2])) {
    if (length(taxon.split.unlist) == 2) return(taxon)
    matches <- unlist(lapply(taxon.split, function(x) x[3:length(x)] %in% words))
    res <- paste(taxon.split.unlist[c(TRUE, TRUE, matches)], collapse = " ")
  } else {
    matches <- unlist(lapply(taxon.split, function(x) x[2:length(x)] %in% words))
    res <- paste(taxon.split.unlist[c(TRUE, matches)], collapse = " ")
  }
  if (length(ident) == 0L) {
    return(res)
  } else {
    return(paste(res, trim(ident)))
  }
}