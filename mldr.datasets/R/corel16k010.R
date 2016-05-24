#'  Datasets with data from the Corel image collection. There are 10 subsets in corel16k
#'
#' @description Multilabel dataset from the image domain.
#' @format An mldr object with 13618 instances, 500 attributes and 144 labels
#' @source Barnard, K. and Duygulu, P. and Forsyth, D. and de Freitas, N. and Blei, D. M. and Jordan, M. I., "Matching words and pictures", Journal of Machine Learning Research, Vol. 3, pp. 1107--1135, 2003
#' @examples
#'\dontrun{
#' corel16k010()  # Check and load the dataset
#' toBibtex(corel16k010)
#' corel16k010$measures
#' }
#' @export
corel16k010 <- function() check_n_load.mldr('corel16k010')
