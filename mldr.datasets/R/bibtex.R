#'  Dataset with BibTeX entries
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 7395 instances, 1836 attributes and 159 labels
#' @source Katakis, I. and Tsoumakas, G. and Vlahavas, I., "Multilabel Text Classification for Automated Tag Suggestion", in Proc. ECML PKDD08 Discovery Challenge, Antwerp, Belgium, pp. 75-83, 2008
#' @examples
#'\dontrun{
#' bibtex()  # Check and load the dataset
#' toBibtex(bibtex)
#' bibtex$measures
#' }
#' @export
bibtex <- function() check_n_load.mldr("bibtex")
