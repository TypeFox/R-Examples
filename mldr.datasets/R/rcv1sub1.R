#'  Dataset from the Reuters corpus (subset 1)
#'
#' @description Multilabel dataset from the text domain.
#' @format An mldr object with 6000 instances, 47236 attributes and 101 labels
#' @source Lewis, D. D. and Yang, Y. and Rose, T. G. and Li, F., "RCV1: A new benchmark collection for text categorization research", The Journal of Machine Learning Research, Vol. 5, pp. 361-397, 2004
#' @examples
#'\dontrun{
#' rcv1sub1()  # Check and load the dataset
#' toBibtex(rcv1sub1)
#' rcv1sub1$measures
#' }
#' @export
rcv1sub1 <- function() check_n_load.mldr('rcv1sub1')
