#' Look up species IDs by taxonomic or common name
#'
#' @export
#' @param name A scientific or common name
#' @param type One of common_name, genus, epithet, or genus_epithet
#' @param fuzzy One of TRUE or FALSE, if FALSE, uses fuzzy search via agrep, if
#'    FALSE, uses grep
#' @examples \dontrun{
#' lookup_names(name='Pinus', type='genus')
#' lookup_names(name='pine', type='common_name')
#' lookup_names(name='bird', type='common_name', fuzzy=TRUE)
#' }
lookup_names <- function(name, type = 'genus', fuzzy = FALSE) {
  type <- match.arg(type, choices = c('common_name','genus','epithet','genus_epithet'))
  if (fuzzy) {
    taxonlist[agrep(name, taxonlist[, type]), ]
  } else {
    taxonlist[grep(name, taxonlist[, type]), ]
  }
}
