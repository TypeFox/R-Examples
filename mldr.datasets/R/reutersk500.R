#' Dataset from the Reuters Corpus with the 500 most relevant features selected
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 6000 instances, 500 attributes and 103 labels
#' @source Read, Jesse, "Scalable multi-label classification", University of Waikato, 2010
#' @examples
#'\dontrun{
#' reutersk500()  # Check and load the dataset
#' toBibtex(reutersk500)
#' reutersk500$measures
#' }
#' @export
reutersk500 <- function() check_n_load.mldr('reutersk500')
