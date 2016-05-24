#'  Dataset from the Stack Exchange's cooking forum
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 10491 instances, 577 attributes and 400 labels
#' @source Charte, Francisco and Rivera, Antonio J. and del Jesus, Maria J. and Herrera, Francisco, "QUINTA: A question tagging assistant to improve the answering ratio in electronic forums", in EUROCON 2015 - International Conference on Computer as a Tool (EUROCON), IEEE, pp. 1-6, 2015
#' @examples
#'\dontrun{
#' stackex_cooking()  # Check and load the dataset
#' toBibtex(stackex_cooking)
#' stackex_cooking$measures
#' }
#' @export
stackex_cooking <- function() check_n_load.mldr('stackex_cooking')
