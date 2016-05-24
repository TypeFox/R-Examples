#'  Dataset from the Stack Exchange's chemistry forum
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 6961 instances, 540 attributes and 175 labels
#' @source Charte, Francisco and Rivera, Antonio J. and del Jesus, Maria J. and Herrera, Francisco, "QUINTA: A question tagging assistant to improve the answering ratio in electronic forums", in EUROCON 2015 - International Conference on Computer as a Tool (EUROCON), IEEE, pp. 1-6, 2015
#' @examples
#'\dontrun{
#' stackex_chemistry()  # Check and load the dataset
#' toBibtex(stackex_chemistry)
#' stackex_chemistry$measures
#' }
#' @export
stackex_chemistry <- function() check_n_load.mldr('stackex_chemistry')
