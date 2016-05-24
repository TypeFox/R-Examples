#'  List  with 10 folds of the train data from the EUR-Lex directory codes dataset
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 17413 instances, 5000 attributes and 412 labels
#' @source Mencia, E. L. and Furnkranz, J., "Efficient pairwise multilabel classification for large-scale problems in the legal domain", Machine Learning and Knowledge Discovery in Databases, pp. 50--65, 2008
#' @examples
#'\dontrun{
#' eurlexdc_tra()  # Check and load the dataset
#' toBibtex(eurlexdc_test[[1]])
#' eurlexdc_test[[1]]$measures
#' }
#' @export
eurlexdc_tra <- function() check_n_load.mldr('eurlexdc_tra')
