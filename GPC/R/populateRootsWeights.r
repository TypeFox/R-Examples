populateRootsWeights <- function(N,L,Z,PolyName,ParamDistrib,NDIi){
  # cat("In popRootsWeights\n")
  
  # dummy= NDIi
  NDXi = rep(0,N*Z)
  dim(NDXi) <- c(N,Z)	
  NDWi = rep(0,Z)
  weighttmp = rep(0,(max(L)*N))
  dim(weighttmp) = c(max(L),N)
  
  for (nn in 1:N){
    if (!is.null(ParamDistrib)){
      if (PolyName[nn] == "JACOBI"){ 
        alpha = ParamDistrib$alpha[nn]
        beta = ParamDistrib$beta[nn]    
      } else if (PolyName[nn] == "LAGUERRE"){
        alpha = ParamDistrib$alpha[nn]
        beta = 0
      } else {stop('Jacobi requires ParamDistrib$alpha and ParamDistrib$beta to be defined.\nLaguerre requires ParamDistrib$shape to be defined.')}
    } else{
      alpha = 0 #rep(0,N)
      beta = 0 #rep(0,N)
    }
    # Computing the roots
    # print(c(L[nn],PolyName[nn],alpha[nn],beta[nn]))
    # coeff = polyTermNorm(L[nn],PolyName[nn],alpha[nn],beta[nn])
    
    # cat(L[nn],PolyName[nn],alpha,beta,'\n')
    coeff = polyTermNorm(L[nn],PolyName[nn],alpha,beta)
    # cat(coeff,'\n')
    
    xi = Re(polyroot(coeff))		
    # Computing the weights
    weight = rep(0.0,L[nn]) 
    for (ii in (1:L[nn])){
      weight[ii] = polyWeight(L[nn],xi,ii,PolyName[nn],alpha,beta)
    }
    weight = weight/sum(weight)
    # Sort the roots and weights with increasing roots		
    ind = sort.list(xi)
    xi = xi[ind]
    weight = weight[ind]		
    NDXi[nn,] = xi[NDIi[nn,]]
    weighttmp[1:L[nn],nn] = weight
  }
  
  for (zz in 1:Z){
    tmp = 1
    for (nn in 1:N){
      tmp = tmp*weighttmp[NDIi[nn,zz],nn]
      NDWi[zz] = tmp
    }
  }
  
  res = list(listNDXi = NDXi, listNDWi = NDWi)
  return(res)
}
