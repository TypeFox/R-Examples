#' Dataset generated from the Yahoo! web site index (reference category)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 8027 instances, 39679 attributes and 33 labels
#' @source Ueda, N. and Saito, K., "Parametric mixture models for multi-labeled text", Advances in neural information processing systems, pp. 721--728, 2002
#' @examples
#'\dontrun{
#' yahoo_reference()  # Check and load the dataset
#' toBibtex(yahoo_reference)
#' yahoo_reference$measures
#' }
#' @export
yahoo_reference <- function() check_n_load.mldr('yahoo_reference')
