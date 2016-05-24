genLargerBAnetwork <-
function(p, sparsity, str=.6){
  ## generate BA networks with over 10 nodes - current SINGLE code is not well suited for this
  #
  #
  #
  # we ensure positive definiteness by removing eigenvectors associated with negative eigenvalues!
  #
  #
  
  
  ## first map the sparsity level to the corresponding value of m:
  sparseLvls = sapply(seq(1,20), FUN=function(m){
    x = as.matrix(get.adjacency(barabasi.game(n = p, power = 1, m = m, directed = FALSE)))
    return(sum(x[upper.tri(x)]!=0)/choose(p,2))
  })
  
  mVal = which.min(abs(sparseLvls-sparsity))
  
  # simulate small world network:
  sub_matrix = as.matrix(get.adjacency(barabasi.game(n = p, power = 1, m = mVal, directed = FALSE)))
  
  # simulate edge weights:
  W = matrix(0, ncol=p, nrow=p)
  weights = runif(p*(p-1)*0.5, str/2, str)
  #W[lower.tri(W)] = weights
  W[upper.tri(W)] = weights
  W = W+t(W)
  sub_matrix = sub_matrix * W
  # simulate edge signs:
  S = matrix(0, ncol=ncol(sub_matrix), nrow=nrow(sub_matrix)) #  matrix of signs
  signs = sample(c(-1,1), p*(p-1)*0.5, replace=TRUE)
  #S[lower.tri(S)] = signs
  S[upper.tri(S)] = signs
  S = S+t(S)
  sub_matrix = sub_matrix * S
  
  # normalize in order to ensure pos definite - we follow Danaher et al 2013
  norm = apply(sub_matrix,1,FUN=function(x){sum(abs(x))})
  sub_matrix = apply(sub_matrix,2, FUN=function(x){x/(1.5*norm)})
  
  sub_matrix = (sub_matrix + t(sub_matrix))/2 # to make symmetric
  diag(sub_matrix) = 1
  
  return(list(TN = sub_matrix, Pres = sub_matrix))
}
