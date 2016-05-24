#'  Dataset with features extracted from video sequences and semantic concepts assigned as labels
#'
#' @description Multilabel dataset from the video domain.
#' @format An mldr object with 43907 instances, 120 attributes and 101 labels
#' @source Snoek, C. G. M. and Worring, M. and van Gemert, J. C. and Geusebroek, J. M. and Smeulders, A. W. M., "The challenge problem for automated detection of 101 semantic concepts in multimedia", in Proc. 14th ACM International Conference on Multimedia, MULTIMEDIA06, pp. 421-430, 2006
#' @examples
#'\dontrun{
#' mediamill()  # Check and load the dataset
#' toBibtex(mediamill)
#' mediamill$measures
#' }
#' @export
mediamill <- function() check_n_load.mldr('mediamill')
