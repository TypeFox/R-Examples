#' Calculate KL divergence of features
#'
#' Computes Kullback-Leibler divergence between features and target vector.
#'
#' @inheritParams calc_ig
#' @return a \code{integer} vector of length equal to the number of features 
#' containing computed information gain values.
#' @note Both \code{target} and \code{features} must be binary, i.e. contain only 0 
#' and 1 values.
#' @seealso \code{\link{test_features}}.
#' 
#' Kullback-Leibler divergence is calculated using \code{\link[entropy]{KL.plugin}}.
#' @export
#' @references STH here
#' @examples tar <- sample(0L:1, 100, replace = TRUE)
#' feat <- sample(0L:1, 100, replace = TRUE)
#' prop <- c(100 - sum(tar), sum(tar))/100
#' entr <- - sum(prop*log(prop))
#' library(bit) #used to code vector as bit
#' calc_kl(feat, as.bit(tar), 100, sum(tar), entr)
calc_kl <- function(feature, target_b, len_target, pos_target, ES) {
  crosstable <- fast_crosstable(target_b, len_target, pos_target, feature)
  counts_feature <- c(crosstable[2] + crosstable[4], crosstable[1] + crosstable[3])
  
  KL.plugin(crosstable[c(2, 4)]/counts_feature[1], crosstable[c(1, 3)]/counts_feature[2])
}