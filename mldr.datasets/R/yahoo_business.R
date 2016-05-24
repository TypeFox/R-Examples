#' Dataset generated from the Yahoo! web site index (business category)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 11214 instances, 21924 attributes and 30 labels
#' @source Ueda, N. and Saito, K., "Parametric mixture models for multi-labeled text", Advances in neural information processing systems, pp. 721--728, 2002
#' @examples
#'\dontrun{
#' yahoo_business()  # Check and load the dataset
#' toBibtex(yahoo_business)
#' yahoo_business$measures
#' }
#' @export
yahoo_business <- function() check_n_load.mldr('yahoo_business')
