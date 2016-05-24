#' Dataset generated from the Yahoo! web site index (society category)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 14512 instances, 31802 attributes and 27 labels
#' @source Ueda, N. and Saito, K., "Parametric mixture models for multi-labeled text", Advances in neural information processing systems, pp. 721--728, 2002
#' @examples
#'\dontrun{
#' yahoo_society()  # Check and load the dataset
#' toBibtex(yahoo_society)
#' yahoo_society$measures
#' }
#' @export
yahoo_society <- function() check_n_load.mldr('yahoo_society')
