#'  Dataset generated from the IMDB film database
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 120919 instances, 1001 attributes and 28 labels
#' @source Read, J. and Pfahringer, B. and Holmes, G. and Frank, E., "Classifier chains for multi-label classification", Machine Learning, (3)85, pp. 333-359, 2011
#' @examples
#'\dontrun{
#' imdb()  # Check and load the dataset
#' toBibtex(imdb)
#' imdb$measures
#' }
#' @export
imdb <- function() check_n_load.mldr('imdb')
