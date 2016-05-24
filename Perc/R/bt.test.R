#' Systemic test for the assumptions of the Bradley-Terry model
#' 
#' \code{bt.test} Systemic test for the assumptions of the Bradley-Terry model,
#'  transitivity and monotonic win-loss relationship. 
#'  That is, if \eqn{A > B} and \eqn{B > C} then \eqn{A > C} and 
#'  \eqn{Pr(A beats C)} > \eqn{Pr(B beats C)}.
#' 
#' @param conf.mat an N-by-N matrix. The matrix should be a conflict matrix with element i,j
#' representing the number of times i has beaten j.
#' @param baseline an integer between 1 and N inclusive identifying 
#' the agent with dominance index equal to zero.
#' @param maxLength an integer indicating maximum path length used
#'  in \code{conductance}
#' @param reps an integer indicating number of conflict matrices
#'  simulated to estimate the sampling distribution under the BT model.
#' @return A list of 3 elements. 
#'  \item{stat}{value of the test statistic}
#'  \item{dist}{estimated sampling distribution of the test statistics under the BT model.}
#'  \item{p.val}{p-value of the test}
#' 
#' @details The value of the test statistic should be within the estimated
#' sampling distribution of the test statistics under the BT model. 
#' The p-value of the test indicates the probability of statistics in the estimated sampling distribution
#' is larger than the test statistic. 
#' It is not appropriate to use Bradley-Terry model if value of the test statistic is higher than 
#' the estimated sampling distribution of the test statistics.
#' 
#' @references 
#'  
#'  Shev, A., Fujii, K., Hsieh, F., & McCowan, B. (2014). Systemic Testing on Bradley-Terry Model against Nonlinear Ranking Hierarchy. PloS one, 9(12), e115367.
#'  
#' @examples
#' # create an edgelist
#' edgelist1 <- data.frame(col1 = sample(letters[1:15], 200, replace = TRUE),
#'                         col2 = sample(letters[1:15], 200, replace = TRUE),
#'                        stringsAsFactors = FALSE)
#' edgelist1 <- edgelist1[-which(edgelist1$col1 == edgelist1$col2), ]
#' # convert an edgelist to conflict matrix
#' confmatrix_bt <- as.conflictmat(edgelist1)
#' # test the assumptions of the Bradley-Terry model
#' # not run:
#' # condTestoutput <- bt.test(confmatrix_bt)
#' @export


###############################################################################
###Description: Systemic test for the assumptions of the Bradley-Terry model,
###             transitivity and monotonic dominance. That is, if A > B and
###             B > C then A > C and Pr(A beats C) > Pr(B beats C).
###Input: 
###  conf.mat - N-by-N conflict matrix whose (i,j)th element is the number of 
###             times i defeated j
###  baseline - integer between 1 and N inclusive identifying the agent with
###             dominance index equal to zero.
###  maxLength - maximum path length used in conductance
###  reps - number of conflict matrices simulated to estimate the sampling
###         distribution under the BT model.
###Ouput: List of 3 items,
###  $stat - value of the test statistic.
###  $dist - estimated sampling distribution of the test statistics under the
###          BT model.
###  $p.val - p-value of the test.
###############################################################################
bt.test = function(conf.mat, baseline = 1, maxLength = 2, reps = 1000){  # temp revision, just to run fastly
  
  # making sure conf is of conf.mat
  if (!("conf.mat" %in% class(conf.mat))){
    conf.mat = as.conflictmat(conf.mat)
  }
  n = nrow(conf.mat)
  mle.d = bradleyTerry(conf.mat, baseline = baseline)
  mle.probs = convertToProb(mle.d[[1]]) # use only domInds (vector of length N consiting of the MLE values of the dominance indices.)
  mle.ord = order(mle.d[[1]], decreasing = TRUE)
  cond = conductance(conf.mat, maxLength)
  test.stat = bt.cond.dist(cond$p.hat, mle.probs, mle.ord)
  num.comps = matrix(0, nrow = n, ncol = n)
  num.comps[upper.tri(num.comps)] = conf.mat[upper.tri(conf.mat)] + 
    t(conf.mat)[upper.tri(conf.mat)]
  num.comps[lower.tri(num.comps)] = t(num.comps)[lower.tri(num.comps)] 
  test.dist = replicate(reps, sampleDist(mle.probs, num.comps, baseline, 
                                         maxLength))  # distance not used
  return(list(stat = test.stat, dist = test.dist, 
              p.val = sum((test.dist>test.stat) / reps)))
}
