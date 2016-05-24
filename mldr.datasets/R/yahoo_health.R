#' Dataset generated from the Yahoo! web site index (health category)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 8205 instances, 30605 attributes and 32 labels
#' @source Ueda, N. and Saito, K., "Parametric mixture models for multi-labeled text", Advances in neural information processing systems, pp. 721--728, 2002
#' @examples
#'\dontrun{
#' yahoo_health()  # Check and load the dataset
#' toBibtex(yahoo_health)
#' yahoo_health$measures
#' }
#' @export
yahoo_health <- function() check_n_load.mldr('yahoo_health')
