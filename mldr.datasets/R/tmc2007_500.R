#'  Dataset from airplanes failures reports (500 most relevant features extracted)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 28596 instances, 500 attributes and 22 labels
#' @source Srivastava, A. N. and Zane-Ulman, B., "Discovering recurring anomalies in text reports regarding complex space systems", Aerospace Conference, pp. 3853-3862, 2005
#' @examples
#'\dontrun{
#' tmc2007_500()  # Check and load the dataset
#' toBibtex(tmc2007_500)
#' tmc2007_500$measures
#' }
#' @export
tmc2007_500 <- function() check_n_load.mldr('tmc2007_500')
