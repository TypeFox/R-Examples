#' Get downstream taxa
#' 
#' Get all downstream taxa from a family or genus name.
#'
#' @param taxon a character vector with either a family or genus name
#' @param accepted list only accepted names?
#' @export
#' @examples
#' \dontrun{
#' lower.taxa("Acosmium")
#' lower.taxa("Zygophyllaceae")
#' }

lower.taxa <- function(taxon, accepted = TRUE) {
  taxon <- fixCase(trim(taxon))
  matches <- apply(all.taxa[, c("family", "genus")], 2, function(x) grepl(paste("^", taxon, "$", sep = ""), x))
  matches <- rowSums(matches) == 1L
  if (accepted) {
    all.taxa$search.str[with(all.taxa, matches & taxon.status == "accepted" & !is.na(taxon.status))]
  } else {
    all.taxa$search.str[matches]
  }
}