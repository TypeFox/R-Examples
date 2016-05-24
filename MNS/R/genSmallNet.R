genSmallNet <-
function(p, pow, m, strUpper, strLower){
  # generate network for one subcomponent as in Danaher et al. 2014
  #
  # INPUT:
  #      - p: number of nodes
  #      - pow, m: power of preferential attachment and number of edges to add at each step (from barabasi.game function in igraph)
  #      - strUpper, strLower: define interval from which to sample edge strengths
  #
  # OUTPUT:
  #      - Adj: adjacency matrix for precision matrix
  #      - Sigma: sampling covariance matrix
  #
  
  posDefInd = FALSE
  
  # get adj matrix:
  while (!posDefInd){
    Adj = as.matrix(get.adjacency(barabasi.game(n = p, power = pow, m = m, directed = FALSE)))
    
    # randomly sample uniform weights:
    U = matrix(0, ncol=p, nrow=p)
    S = matrix(0, ncol=p, nrow=p)
    
    U[upper.tri(U)] = runif(choose(p,2), min = strLower, max = strUpper)
    U = U + t(U)
    S[upper.tri(S)] = sample(c(-1,1), size=choose(p,2), replace=TRUE)
    S = S + t(S)
    
    #U = matrix(runif(p*p, min = strLower, max = strUpper), ncol=p)
    #S = matrix(sample(c(-1,1), size = p*p, replace = TRUE), ncol=p)
    #S = matrix(0, ncol=p, nrow=p)
    #S[upper.tri(S)] = sample(c(-1,1), size=choose(p,2), replace = TRUE)
    #S = S+t(S)
    
    Adj = Adj* U * S
    
    # normalize in order to ensure pos definite:
    norm = apply(Adj,1,FUN=function(x){sum(abs(x))})
    Adj = apply(Adj,2, FUN=function(x){x/(1.5*norm)})
    
    Adj = (Adj + t(Adj))/2 # to make symmetric
    diag(Adj) = 1
    
    # now get sampling covariance:
    Sigma = solve(Adj)
    diag(Sigma) = abs(diag(Sigma))
    Sigma = Sigma/sqrt(abs(diag(Sigma)))%*%sqrt(t(abs(diag(Sigma))))
    
    Sigma = Sigma * .6
    diag(Sigma) = 1
    
    # check we can generate MVG to ensure its pos def!
    s = try(mvrnorm(n=1, mu = rep(0,p), Sigma = Sigma), silent = TRUE)
    if (class(s)!="try-error") posDefInd=TRUE
  }
  
  return(list(Adj=Adj, Sigma=Sigma))
  
}
