#'  Dataset from airplanes failures reports
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 28596 instances, 49060 attributes and 22 labels
#' @source Srivastava, A. N. and Zane-Ulman, B., "Discovering recurring anomalies in text reports regarding complex space systems", Aerospace Conference, pp. 3853-3862, 2005
#' @examples
#'\dontrun{
#' tmc2007()  # Check and load the dataset
#' toBibtex(tmc2007)
#' tmc2007$measures
#' }
#' @export
tmc2007 <- function() check_n_load.mldr('tmc2007')
