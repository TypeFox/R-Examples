#' Dataset generated from the Yahoo! web site index (science category)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 6428 instances, 37187 attributes and 40 labels
#' @source Ueda, N. and Saito, K., "Parametric mixture models for multi-labeled text", Advances in neural information processing systems, pp. 721--728, 2002
#' @examples
#'\dontrun{
#' yahoo_science()  # Check and load the dataset
#' toBibtex(yahoo_science)
#' yahoo_science$measures
#' }
#' @export
yahoo_science <- function() check_n_load.mldr('yahoo_science')
