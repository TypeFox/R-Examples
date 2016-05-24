#'  Dataset with data from the Corel image collection
#'
#' @description Multilabel dataset from the image domain.
#' @format An mldr object with 5000 instances, 499 attributes and 374 labels
#' @source Duygulu, P. and Barnard, K. and de Freitas, J.F.G. and Forsyth, D.A., "Object Recognition as Machine Translation: Learning a Lexicon for a Fixed Image Vocabulary", Computer Vision, ECCV 2002, pp. 97-112, 2002
#' @examples
#'\dontrun{
#' corel5k()  # Check and load the dataset
#' toBibtex(corel5k)
#' corel5k$measures
#' }
#' @export
corel5k <- function() check_n_load.mldr("corel5k")
