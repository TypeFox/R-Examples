NetGSA <-
function(
  A1,     		  
  A2, 			    	
  x, 			      	
  y, 			      	
  B, 		    	  	
  lklMethod = c("REML","ML"),
  directed = FALSE,          
  eta = 1e-1,           
  lim4kappa = 500       
){
  this.call <- match.call()
  lklMethod <- match.arg(lklMethod)
  # s2profile <- match.arg(s2profile)
  
  p = dim(x)[1] #No. of genes
  n = length(y) #No. of samples in total
  
  if (dim(x)[2] != n) {
    stop("The dimensions of the data matrix and class vector don't match.")
  }
  
  if (dim(B)[2] != p) {
    stop("The dimensions of the data matrix and indicator matrix don't match.")
  }
  
  if (length(unique(y)) != 2) {
    stop("There should be 2 unique classes in the class indicator!")
  }
  
  if ((sum(abs(A1)) == 0) && (sum(abs(A2)) == 0)) {
    stop("No network interactions were found!")
  }
  
  ##-----------------
  ##Determine whether the network is DAG
  ##Assume A1 and A2 are of the same type (directed or undirected)
  isNetDAG = FALSE
  if (directed) {    
	  gA <- graph.adjacency(A1, mode="directed")
	  isNetDAG = is.dag(gA)
  }  

  if (p > 5000) {
    warning("netGSA may be slow for datasets with large number of genes.")
  }
  
  ##-----------------
  ##setting up control parameters for the var estimation procedures
  varEstCntrl = list(lklMethod = lklMethod,                    
                     s2profile = "se",   
                     lb = 0.5,           
                     ub = 100,           
                     tol = 0.01)         
  
  ##-----------------
  ##Find influence matrices based on adjacency matrices A1 and A2
  if (directed){
    D1 = adj2inf(AA = A1, isDAG = isNetDAG, eta = eta)
    D2 = adj2inf(AA = A2, isDAG = isNetDAG, eta = eta)
    
    tmp = min(kappa(D1), kappa(D2))
    
    while ((min(kappa(D1), kappa(D2)) > lim4kappa) && !isNetDAG) {
      eta = eta * 2
      D1 = adj2inf(AA = A1, isDAG = isNetDAG, eta = eta)
      D2 = adj2inf(AA = A2, isDAG = isNetDAG, eta = eta)
    }
    
    if (((tmp > lim4kappa) && !isNetDAG)) {
      warning(paste("Influence matrix is ill-conditioned, using eta =", eta))
    } 
    
    D1D1 = D1 %*% t(D1)
    D2D2 = D2 %*% t(D2)
    
    tmp = min(kappa(D1D1), kappa(D2D2))
    while ((tmp > lim4kappa) && !isNetDAG) {
      eta = eta * 2
      
      D1 = adj2inf(AA = A1, isDAG = isNetDAG, eta = eta)
      D2 = adj2inf(AA = A2, isDAG = isNetDAG, eta = eta)
      D1D1 = D1 %*% t(D1)
      D2D2 = D2 %*% t(D2)
    }
    
    if (((tmp > lim4kappa) && !isNetDAG)) {
      warning(paste("Influence matrix is ill-conditioned, using eta =", eta))  
    }
  } else {
    mat4norm = normalizeAdj(list(A1, A2), alpha=0.99)
    
    D1 = mat4norm$InfMat[[1]]
    D2 = mat4norm$InfMat[[2]]
  }
  
  output = call.netGSA(D1, D2, x, y, B, varEstCntrl)
  
  return(output)
}
