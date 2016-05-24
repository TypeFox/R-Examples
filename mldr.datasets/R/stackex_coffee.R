#'  Dataset from the Stack Exchange's coffee forum
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 225 instances, 1763 attributes and 123 labels
#' @source Charte, Francisco and Rivera, Antonio J. and del Jesus, Maria J. and Herrera, Francisco, "QUINTA: A question tagging assistant to improve the answering ratio in electronic forums", in EUROCON 2015 - International Conference on Computer as a Tool (EUROCON), IEEE, pp. 1-6, 2015
#' @examples
#'\dontrun{
#' stackex_coffee()
#' toBibtex(stackex_coffee)
#' stackex_coffee$measures
#' }
#' @export
stackex_coffee <- function() check_n_load.mldr('stackex_coffee')
