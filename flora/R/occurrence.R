#' Taxa occurrence
#' 
#' Find the taxa that occur in a given state of Brazil.
#' 
#' @param states a character vector with one or more state abbreviations
#'   following. See notes for abbreviations.
#' @param type type of matching to be used. \code{any} will return the taxa that
#'   occur in any of the passed \code{states}. \code{only} matches taxa that
#'   occur only in all provided (no more, no less) \code{states} and \code{all} matches taxa that
#'   occur at least in all \code{states} passed. See examples.
#' @param taxa optional character vector to match against the states
#' @export
#' @note List of abbreviations: \url{http://en.wikipedia.org/wiki/States_of_Brazil}
#' @return a data frame
#' @examples
#' \dontrun{
#' occ.any <- occurrence(c("SP", "BA", "MG"), type = "any")
#' occ.only <- occurrence(c("SP", "BA", "MG"), type = "only")
#' occ.all <- occurrence(c("SP", "BA", "MG"), type = "all")
#' occ.taxa <- occurrence(c("SP", "BA", "MG"), type = "all", taxa = lower.taxa("Myrcia"))
#' 
#' head(occ.any)
#' head(occ.only)
#' head(occ.all)
#' head(occ.taxa)
#' }
occurrence <- function(states, type = c("any", "only", "all"), taxa = NULL) {
  type <- match.arg(type)
  states <- sort(sapply(trim(states), toupper))
  #res <- lapply(occurrences, match, states)
  if (type == "any") {
    #res <- lapply(res, function(x) any(!is.na(x)))
    res <- subset(distribution, grepl(paste(states, collapse = "|"), occurrence))
  }
  if (type == "only") {
    res <- subset(distribution, grepl(paste("^", paste(states, collapse = ";"), "$", sep = ""), occurrence))
  }
  if (type == "all") {
    res <- subset(distribution, grepl(paste(states, collapse = ".*"), occurrence))
  }
  # res <- distribution[unlist(res), ]
  if (nrow(res) == 0) {
    return(NA)
  }
  if (is.null(taxa)) {
    merge(all.taxa[, c("id", "family", "search.str")], res[, c("id", "occurrence")], by = "id")
  } else {
    merge(all.taxa[all.taxa$search.str %in% taxa, c("id", "family", "search.str")], res[, c("id", "occurrence")], by = "id")
  }
}