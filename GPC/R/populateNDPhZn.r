populateNDPhZn <- function(N,M,Z,PolyName,NDXi,Index,ParamDistrib){
  # alpha = ParamDistrib$alpha
  # beta = ParamDistrib$beta
  NDPhZn = rep(0,M*Z)
  dim(NDPhZn) <- c(M,Z)	
  for (zz in 1:Z){
    for (mm in 1:M){
      tmp = 1
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
        # tmp = tmp*polyNorm(Index[nn,mm],NDXi[nn,zz],PolyName[nn],alpha[nn],beta[nn])
        tmp = tmp*polyNorm(Index[nn,mm],NDXi[nn,zz],PolyName[nn],alpha,beta)
        NDPhZn[mm,zz] = tmp    	
      }			
    }
  }
  # print(NDPhZn)
  return(NDPhZn)
}