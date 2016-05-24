#' Dataset generated from the Yahoo! web site index (arts entertainment)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 12730 instances, 32001 attributes and 21 labels
#' @source Ueda, N. and Saito, K., "Parametric mixture models for multi-labeled text", Advances in neural information processing systems, pp. 721--728, 2002
#' @examples
#'\dontrun{
#' yahoo_entertainment()  # Check and load the dataset
#' toBibtex(yahoo_entertainment)
#' yahoo_entertainment$measures
#' }
#' @export
yahoo_entertainment <- function() check_n_load.mldr('yahoo_entertainment')
