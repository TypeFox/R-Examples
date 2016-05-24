#' Dataset generated from the Yahoo! web site index (arts category)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 7484 instances, 23146 attributes and 26 labels
#' @source Ueda, N. and Saito, K., "Parametric mixture models for multi-labeled text", Advances in neural information processing systems, pp. 721--728, 2002
#' @examples
#'\dontrun{
#' yahoo_arts()  # Check and load the dataset
#' toBibtex(yahoo_arts)
#' yahoo_arts$measures
#' }
#' @export
yahoo_arts <- function() check_n_load.mldr('yahoo_arts')
