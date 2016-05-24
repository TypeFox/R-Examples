#' Dataset generated from the Yahoo! web site index (recreation category)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 12828 instances, 30324 attributes and 22 labels
#' @source Ueda, N. and Saito, K., "Parametric mixture models for multi-labeled text", Advances in neural information processing systems, pp. 721--728, 2002
#' @examples
#'\dontrun{
#' yahoo_recreation()  # Check and load the dataset
#' toBibtex(yahoo_recreation)
#' yahoo_recreation$measures
#' }
#' @export
yahoo_recreation <- function() check_n_load.mldr('yahoo_recreation')
