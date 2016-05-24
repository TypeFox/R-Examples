#' Dataset generated from the Yahoo! web site index (social category)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 12111 instances, 52350 attributes and 39 labels
#' @source Ueda, N. and Saito, K., "Parametric mixture models for multi-labeled text", Advances in neural information processing systems, pp. 721--728, 2002
#' @examples
#'\dontrun{
#' yahoo_social()  # Check and load the dataset
#' toBibtex(yahoo_social)
#' yahoo_social$measures
#' }
#' @export
yahoo_social <- function() check_n_load.mldr('yahoo_social')
