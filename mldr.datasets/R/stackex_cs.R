#'  Dataset from the Stack Exchange's computer science forum
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 9270 instances, 635 attributes and 274 labels
#' @source Charte, Francisco and Rivera, Antonio J. and del Jesus, Maria J. and Herrera, Francisco, "QUINTA: A question tagging assistant to improve the answering ratio in electronic forums", in EUROCON 2015 - International Conference on Computer as a Tool (EUROCON), IEEE, pp. 1-6, 2015
#' @examples
#'\dontrun{
#' stackex_cs()  # Check and load the dataset
#' toBibtex(stackex_cs)
#' stackex_cs$measures
#' }
#' @export
stackex_cs <- function() check_n_load.mldr('stackex_cs')
