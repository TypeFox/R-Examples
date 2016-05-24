#'  Datasets with data from the Corel image collection. There are 10 subsets in corel16k
#'
#' @description Multilabel dataset from the image domain.
#' @format An mldr object with 13760 instances, 500 attributes and 154 labels
#' @source Barnard, K. and Duygulu, P. and Forsyth, D. and de Freitas, N. and Blei, D. M. and Jordan, M. I., "Matching words and pictures", Journal of Machine Learning Research, Vol. 3, pp. 1107--1135, 2003
#' @examples
#'\dontrun{
#' corel16k003()  # Check and load the dataset
#' toBibtex(corel16k003)
#' corel16k003$measures
#' }
#' @export
corel16k003 <- function() check_n_load.mldr('corel16k003')
