#'  List  with 10 folds of the train data from the EUR-Lex EUROVOC descriptors dataset
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 17413 instances, 5000 attributes and 3993 labels
#' @source Mencia, E. L. and Furnkranz, J., "Efficient pairwise multilabel classification for large-scale problems in the legal domain", Machine Learning and Knowledge Discovery in Databases, pp. 50--65, 2008
#' @examples
#'\dontrun{
#' eurlexev_tra()  # Check and load the dataset
#' toBibtex(eurlexev_tra[[1]])
#' eurlexev_tra[[1]]$measures
#' }
#' @export
eurlexev_tra <- function() check_n_load.mldr('eurlexev_tra')
