
###############################################################################
###Description: Computes the log-likelihood for the BT model.
###Input:
###  conf.mat - N-by-N conflict matrix whose (i,j)th element is the number of 
###             times i defeated j
###  d - numeric vector of length N consisting of dominance indices
###  sumToOne - logical, TRUE indicate 'sum-to-one' parameterization of BT
###Output: numeric value representing the log-likelihood for given at d
###############################################################################
BTLogLik = function(conf.mat, d, sumToOne = FALSE){
  n = nrow(conf.mat)
  m1 = matrix(rep(d, times = n), nrow = n)
  m2 = matrix(rep(d, times = n), nrow = n, byrow = TRUE)
  log.lik = sum(conf.mat * (m1 - log(exp(m1) + exp(m2))))
  return(log.lik)
}




###############################################################################
###Description: Sample a value of the "systemic test" test statistic for data
###             generated under the BT model.
###Input:
###  prob.mle - N-by-N numeric matrix who entry [i,j] is the BT estimate for
###             the probability that i is dominant over j in a game.
###  num.comps - N-by-N numeric matrix whose entry, [i,j], is the number
###              of observed games between agent i and agent j.
###  baseline - interger between 1 and N inclusive.  ID of baseline agent with
###             dominance index set to 0.
###  maxLength - integer, maximum length of paths used in conductance.
###Output: numeric, value of test statistic for data sampled under BT
###############################################################################
sampleDist = function(prob.mle, num.comps, baseline, maxLength){
  n = nrow(prob.mle)  
  conf = sampleBTConf(n, prob.mle, num.comps) # varcov not used
  conf.bt = bradleyTerry(conf, baseline = baseline)
  conf.ord = order(conf.bt[[1]], decreasing = TRUE)  # conf.bt --> conf.bt[[1]] - Error in order(conf.bt, decreasing = TRUE) : unimplemented type 'list' in 'orderVector1'
  conf.cond = conductance(conf, maxLength)
  d = bt.cond.dist(conf.cond$p.hat, convertToProb(conf.bt[[1]]), conf.ord) ## conf.bt --> conf.bt[[1]]
  return(d)
}


###############################################################################
###Description: Sample a conflict matrix based on probabilities from a BT 
###             model.
###Input:
###  n - number of agents (N)
###  p.mat - N-by-N numeric matrix who entry [i,j] is the BT estimate for
###             the probability that i is dominant over j in a game.
###  num.comps - N-by-N numeric matrix whose entry, [i,j], is the number
###              of observed games between agent i and agent j.
###Output: N-by-N conflict matrix simulated under the BT model
###############################################################################
sampleBTConf = function(n, p.mat, num.comps){
  conf = matrix(0, n, n)
  conf[upper.tri(conf)] = rbinom(.5 * n * (n - 1), 
                                 num.comps[upper.tri(num.comps)],
                                 p.mat[upper.tri(p.mat)])
  conf[lower.tri(conf)] = t(num.comps)[lower.tri(num.comps)] - 
    t(conf)[lower.tri(conf)]
  #if((sum(rowSums(conf) == 0) + sum(colSums(conf) == 0)) > 0){
  #  sampleBTConf(n, p.mat, num.comps)
  #}
  return(conf)
}


###############################################################################
###Description: Computes the distance between two matrices.  Here, the distance
###             is defined as the sum of the square root of the absolute value
###             of the components of each matrix multiplied by the square root
###             of the absolute value of the difference in the row and column
###             number.  The idea is if the matrices are ordered by the ranks
###             of the individuals a bigger differences between high and low
###             ranked agents should have a larger effect.
###Input:
###  m1 - N-by-N numeric matrix of dominance probabilities.
###  m2 - N-by_N numeric matrix of dominance probabilities.
###  ord - integer vector containing integers 1 through N.  Denotes the way
###        the rows and columns of m1 and m2 should be ordered.  This might
###        be the outpout of order(d) where d is dominance indices from the
###        BT model.
###Output: numeric value of the distance between the matrices.
###############################################################################
bt.cond.dist = function(m1, m2, ord){
  mat1 = m1[ord,ord]
  mat2 = m2[ord,ord]
  weights = row(mat1) - col(mat1)
  cost = sum(sqrt(abs(mat1 - mat2)) * sqrt(abs(weights)))
  return(cost)
}



###############################################################################
### Description: function converts BT dominance indices into dominance a 
###              matrix of dominance probabilities
### Input: 
###   d - numeric vector of dominance probabilities estimated by the BT model.
###   sumToOne- logical, if TRUE indices are for the "sum to one" 
###               parameterization of the BT model
### Output: N-by-N numeric matrix of dominance probabilities estimated by the
###         BT model
###############################################################################

convertToProb = function(d, sumToOne = FALSE){
  n = length(d)
  m1 = matrix(rep(d, times = n), nrow = n)
  m2 = matrix(rep(d, times = n), nrow = n, byrow = TRUE)
  if(sumToOne){
    P = m1 / (m1 + m2)
  }
  else{
    P = exp(m1) / (exp(m1) + exp(m2))
  }
  diag(P) = 0
  return(P)
}
