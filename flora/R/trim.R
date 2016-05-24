#' Trim a name and remove duplicate tabs and whitespaces
#' 
#' Remove duplicate and misplaced whitespace characters
#' 
#' @param taxon a character vector with a single taxon name
#' @export
#' @return a character vector
#' @examples
#' \dontrun{
#' trim("   Myrcia  lingua")
#' }
trim <- function(taxon) {
  taxon <- gsub("\\t", " ", taxon)
  taxon <- gsub("\\s+", " ", taxon)
  taxon <- gsub("^\\s+", "", taxon)
  taxon <- gsub("\\s+$", "", taxon)
  return(taxon)
}