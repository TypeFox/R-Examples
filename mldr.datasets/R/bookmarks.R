#'  Dataset with data from web bookmarks and their categories
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 87856 instances, 2150 attributes and 208 labels
#' @source Katakis, I. and Tsoumakas, G. and Vlahavas, I., "Multilabel Text Classification for Automated Tag Suggestion", in Proc. ECML PKDD08 Discovery Challenge, Antwerp, Belgium, pp. 75-83, 2008
#' @examples
#'\dontrun{
#' bookmarks()  # Check and load the dataset
#' toBibtex(bookmarks)
#' bookmarks$measures
#' }
#' @export
bookmarks <- function() check_n_load.mldr("bookmarks")
