#' Calculate IG for single feature
#'
#' Computes information gain of single feature and target vector.
#'
#' @param feature feature vector.
#' @param target_b target in bits (as per \code{\link[bit]{as.bit}}).
#' @param len_target length of target vector.
#' @param pos_target number of positive cases in target vector.
#' @param ES numeric value of target entropy.
#' @return a single numeric value - information gain in nats.
#' @details The information gain term is used here (improperly) as a synonym of mutual information. 
#' It is defined as:
#' \deqn{IG(X; Y) = \sum_{y \in Y} \sum_{x \in X} p(x, y) \log \left(\frac{p(x, y)}{p(x) p(y)}  \right)}
#' 
#' In biogram package information gain is calculated using following relationship: 
#' \eqn{IG = E(S) - E(S|F)}
#' @note During calculations \eqn{0 \log 0  = 0}. For justification see References. 
#' 
#' Input looks strange, but the function was designed to be as fast 
#' as possible subroutine of \code{\link{calc_ig}} and generally should not be directly
#' called by user.
#' @references Cover TM, Thomas JA \emph{Elements of Information Theory, 2nd Edition}
#' Wiley, 2006.
#' @export
#' @examples tar <- sample(0L:1, 100, replace = TRUE)
#' feat <- sample(0L:1, 100, replace = TRUE)
#' prop <- c(100 - sum(tar), sum(tar))/100
#' entr <- - sum(prop*log(prop))
#' library(bit) #used to code vector as bit
#' calc_ig(feat, as.bit(tar), 100, sum(tar), entr)
calc_ig <- function(feature, target_b, len_target, pos_target, ES) {
  crosstable <- fast_crosstable(target_b, len_target, pos_target, feature)
  counts_feature <- c(crosstable[2] + crosstable[4], crosstable[1] + crosstable[3])
  
  log_crosstable <- c(entlog(crosstable[1] %/e% counts_feature[2]),
                      entlog(crosstable[3] %/e% counts_feature[2]),
                      entlog(crosstable[2] %/e% counts_feature[1]),
                      entlog(crosstable[4] %/e% counts_feature[1]))
  
  #entropy - conditional entrophy
  ES + (crosstable[1] * log_crosstable[1] +
          crosstable[3] * log_crosstable[2] +
          crosstable[2] * log_crosstable[3] + 
          crosstable[4] * log_crosstable[4])/len_target
}

#logarithm safe for entropy calculation
entlog <- function(x, ...)
  ifelse(x == 0, 0, log(x, ...))

#division safe for entropy calculation
"%/e%" <- function(x, y)
  ifelse(x == 0 && y == 0, 0, x/y)