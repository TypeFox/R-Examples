#' Vernacular name search
#' 
#' Search for taxa using vernacular names
#' 
#' @param name a vernacular name
#' @param exact approximate or exact match?
#' @export
#' @return a data frame of results or NA
#' @examples
#' \dontrun{
#' vernacular("pimenta", exact = TRUE)
#' vernacular("pimenta", exact = FALSE)
#' }

vernacular <- function(name, exact = FALSE) {
  name <- trim(name)
  if (exact) {
    res <- vernacular.names[grep(paste("^", name, "$", sep = ""), vernacular.names$vernacular.name, ignore.case = TRUE), c("id", "vernacular.name", "locality")]
  } else {
    res <- vernacular.names[agrep(name, vernacular.names$vernacular.name, ignore.case = TRUE), c("id", "locality", "vernacular.name")]
  }
  if (nrow(res) == 0L) {
    NA
  } else {
    merge(all.taxa[, c("id", "search.str", "family")], res, by = "id")
  }
}