#'  Dataset from images with different natural scenes
#'
#' @description Multilabel dataset from the image domain.
#' @format An mldr object with 2407 instances, 294 attributes and 6 labels
#' @source Boutell, M. and Luo, J. and Shen, X. and Brown, C., "Learning multi-label scene classification", Pattern Recognition, (9)37, pp. 1757--1771, 2004
#' @examples
#'\dontrun{
#' scene()
#' toBibtex(scene)
#' scene$measures
#' }
scene <- function() check_n_load.mldr('scene')
