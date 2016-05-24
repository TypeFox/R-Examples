fem.sparse <- function(Y,K,maxit,eps,Tinit,model,method='reg',l1,nbit,l2){
  colnames = colnames(Y)
  if (length(l1)!=1 | l1>1){cat('\n','The l1 penalty term is a single figure comprises between 0 and 1','\n'); break}
  # Initialization
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  d = min((K-1),(p-1))

  # New objects
  Lobs = rep(c(-Inf),1,(maxit+1))
  # 	
  # Initialization of T
  T         = Tinit
  V         = fstep.sparse(Y,T,l1,nbit,l2)
  prms      = fem.mstep(Y,V,T,model=model,method=method)
  res.estep = fem.estep(prms,Y,V)
  T         = res.estep$T
  Lobs[1]   = res.estep$loglik
  
  # Main loop
  Linf_new  = Lobs[1]
      for (i in 1:maxit){
      # The three main steps F, M and E
		      V         = fstep.sparse(Y,T,l1,nbit,l2)
		      prms      = fem.mstep(Y,V,T,model=model,method=method)
		      res.estep = fem.estep(prms,Y,V)
		      T         = res.estep$T
		      Lobs[i+1] = res.estep$loglik
    
    # Stop criterion
    if (i>=2){
      acc      = (Lobs[i+1] - Lobs[i]) / (Lobs[i] - Lobs[i-1])
      Linf_old = Linf_new
      Linf_new <- try( Lobs[i] + 1/(1-acc) * (Lobs[i+1] - Lobs[i]))
      if (abs(Linf_new - Linf_old) < eps) {break}
    }
    
  }
  
  # Returning the results
  cls  = max.col(T)
  crit = fem.criteria(Lobs[(i+1)],T,prms,n)
  rownames(V) = colnames
  colnames(V) = paste('U',1:d,sep='')
  res  = list(K=K,cls=cls,P=T,U=V,aic=crit$aic,mean=prms$mean,my=prms$my,prop=prms$prop,D=prms$D,model=prms$model,bic=crit$bic,icl=crit$icl,loglik=Lobs[2:(i+1)],ll=Lobs[i+1],method=method)
  
  res
}

