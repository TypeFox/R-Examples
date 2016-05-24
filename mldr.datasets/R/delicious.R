#'  Dataset generated from the del.icio.us site bookmarks
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 16105 instances, 500 attributes and 983 labels
#' @source Tsoumakas, G. and Katakis, I. and Vlahavas, I., "Effective and Efficient Multilabel Classification in Domains with Large Number of Labels", in Proc. ECML/PKDD Workshop on Mining Multidimensional Data, Antwerp, Belgium, MMD08, pp. 30--44, 2008
#' @examples
#'\dontrun{
#' delicious()  # Check and load the dataset
#' toBibtex(delicious)
#' delicious$measures
#' }
#' @export
delicious <- function() check_n_load.mldr('delicious')
