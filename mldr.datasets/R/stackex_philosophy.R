#'  Dataset from the Stack Exchange's philosophy forum
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 3971 instances, 842 attributes and 233 labels
#' @source Charte, Francisco and Rivera, Antonio J. and del Jesus, Maria J. and Herrera, Francisco, "QUINTA: A question tagging assistant to improve the answering ratio in electronic forums", in EUROCON 2015 - International Conference on Computer as a Tool (EUROCON), IEEE, pp. 1-6, 2015
#' @examples
#'\dontrun{
#' stackex_philosophy()  # Check and load the dataset
#' toBibtex(stackex_philosophy)
#' stackex_philosophy$measures
#' }
#' @export
stackex_philosophy <- function() check_n_load.mldr('stackex_philosophy')
