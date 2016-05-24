#' Computes the MLE for the BT model using an MM algorithm
#' 
#' \code{bradleyTerry} Computes the MLE for the BT model using an MM algorithm
#' 
#' @param conf.mat a matrix of conf.mat class. An N-by-N conflict matrix whose \code{(i,j)}th element is the number of times i defeated j.
#' @param initial initial values of dominance indices for the MM algorithm, if not supplied, the 0 vector will be the inital value.
#' @param baseline index for agent to represent baseline dominance index set to 0.  If NA, the "sum-to-one" parameterization will be used.
#' @param stop.dif numeric value for difference in log likelihood value between iterations.  Used as the convergence criterion for the algorithm.
#' @return A list of length 3. 
#'  \item{domInds}{a vector of length N consiting of the MLE values of the dominance indices. Lower values represent lower ranks.}
#'  \item{probMat}{an N-by-N numeric matrix of win-loss probabilities estimated by the BT model.}
#'  \item{logLik}{the model fit.}
#'
#' @details In order to meet Bradley-Terry assumption, each ID in \code{conf.mat} should have at least one win AND one loss.
#' \code{bradleyTerry} will return an error if no more than one win or loss was found.
#'  
#'  @references 
#'    Shev, A., Hsieh, F., Beisner, B., & McCowan, B. (2012). Using Markov chain Monte Carlo (MCMC) to visualize and test the linearity assumption of the Bradley-Terry class of models. Animal behaviour, 84(6), 1523-1531.
#'    
#'    Shev, A., Fujii, K., Hsieh, F., & McCowan, B. (2014). Systemic Testing on Bradley-Terry Model against Nonlinear Ranking Hierarchy. PloS one, 9(12), e115367.
#'  
#' @examples
#' # create an edgelist
#' edgelist1 <- data.frame(col1 = sample(letters[1:15], 200, replace = TRUE),
#'                         col2 = sample(letters[1:15], 200, replace = TRUE),
#'                        stringsAsFactors = FALSE)
#' edgelist1 <- edgelist1[-which(edgelist1$col1 == edgelist1$col2), ]
#' # convert an edgelist to conflict matrix
#' confmatrix_bt <- as.conflictmat(edgelist1)
#' # Computes the MLE for the BT model
#' bt <- bradleyTerry(confmatrix_bt)
#' @export

###############################################################################
###Description: Computes the MLE for the BT model using an MM algorithm
###Input:
###  conf.mat - N-by-N conflict matrix whose (i,j)th element is the number of 
###             times i defeated j
###  initial - initial values of dominance indices for the MM algorithm, if 
###            not supplied, the 0 vector will be the inital value.
###  baseline - index for agent to represent baseline dominance index set to
###             0.  If NA, the "sum-to-one" parameterization will be used.
###  stop.dif - numeric value for difference in log likelihood value between
###             iterations.  Used as the convergence criterion for the 
###             algorithm.
###Output: 
###  domInds - vector of length N consiting of the MLE values of the dominance
###        indices.
###
###  probMat - N-by-N numeric matrix of dominance probabilities estimated by the
###         BT model
###############################################################################
bradleyTerry = function(conf.mat, initial = NA, baseline = NA, 
                        stop.dif = .001){
  
  #Check row and columns for players who never lost or never won.
  if((sum(rowSums(conf.mat) == 0) > 0) | (sum(colSums(conf.mat) == 0) > 0))  {
    stop("Conflict Matrix does not meet Bradley-Terry assumption.  MLE does not exist.")
  }
  # making sure conf is of conf.mat
  if (!("conf.mat" %in% class(conf.mat))){
    conf.mat = as.conflictmat(conf.mat)
  }

                        
  m = nrow(conf.mat)
  if(length(initial) == 1){
    initial = rep(1/m, times = m)
  }
  d = initial
  W = rowSums(conf.mat)
  N = matrix(0, nrow = m, ncol = m)
  lik.old = BTLogLik(conf.mat, d, TRUE)
  dif = 1
  its = 0
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      N[i,j] = conf.mat[i,j] + conf.mat[j,i]
      N[j,i] = N[i,j]
    }
  }
  while(dif > stop.dif){
    for(i in 1:m){
      C = sapply((1:m)[-i], FUN = function(j, i, N){N[i,j] / (d[i] + d[j])},
                 i = i, N = N)
      d[i] = W[i] * (1 / sum(C))
      d = d/sum(d)
    }
    d[d == 0] = 10^(-43)
    d[d == 1] = 1 - 1^(-43)
    lik.new = BTLogLik(conf.mat, d, TRUE)
    dif = lik.new - lik.old
    lik.old = lik.new
    its = its + 1
  }
  if(!is.na(baseline)){
    d = log(d) - log(d[baseline])
  }
  logLik = lik.new
  attr(logLik, "iterations") = its
  probMat = convertToProb(d, is.na(baseline))
  rownames(probMat) <- rownames(conf.mat)
  colnames(probMat) <- colnames(conf.mat)
  return(list(domInds = d, probMat = probMat, 
              logLik = logLik))  
}

