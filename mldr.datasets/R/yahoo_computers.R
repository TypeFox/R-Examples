#' Dataset generated from the Yahoo! web site index (computers category)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 12444 instances, 34096 attributes and 33 labels
#' @source Ueda, N. and Saito, K., "Parametric mixture models for multi-labeled text", Advances in neural information processing systems, pp. 721--728, 2002
#' @examples
#'\dontrun{
#' yahoo_computers()  # Check and load the dataset
#' toBibtex(yahoo_computers)
#' yahoo_computers$measures
#' }
#' @export
yahoo_computers <- function() check_n_load.mldr('yahoo_computers')
