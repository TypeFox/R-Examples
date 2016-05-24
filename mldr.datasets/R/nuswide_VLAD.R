#'  Dataset obtained from the NUS-WIDE database with cVLAD+ representation
#'
#' @description Multilabel dataset from the image domain.
#' @format An mldr object with 269648 instances, 129 attributes and 81 labels
#' @source Chua, Tat-Seng and Tang, Jinhui and Hong, Richang and Li, Haojie and Luo, Zhiping and Zheng, Yantao, "NUS-WIDE: a real-world web image database from National University of Singapore", in Proc. of the ACM international conference on image and video retrieval, pp. 48, 2009
#' @examples
#'\dontrun{
#' nuswide_VLAD()  # Check and load the dataset
#' toBibtex(nuswide_VLAD)
#' nuswide_VLAD$measures
#' }
#' @export
nuswide_VLAD <- function() check_n_load.mldr('nuswide_VLAD')
