fem.main <-
function(Y,K,init,nstart,maxit,eps,Tinit,model,kernel,method){
  # Initialization
  colnames = colnames(Y)
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  d = min((K-1),(p-1))
  
  # New objects
  Lobs = rep(c(-Inf),1,(maxit+1))
  
  # Initialization of T
  if (init=='user'){ T = Tinit}
  else if (init=='hclust'){
	T   = matrix(0,n,K)
	ind = cutree(hclust(dist(Y),method='ward'),K)
	for (i in 1:n){ T[i,ind[i]] = 1 }
  }
  else if (init=='kmeans' || init=='random'){
	Ltmp = rep(NA,nstart); TT = list()
	for (i in 1:nstart){
		if (init=='random'){TT[[i]] = t(rmultinom(n,1,c(rep(1/K,K))))}
		else{
   			T = matrix(0,n,K)
    			ind = kmeans(Y,K,nstart=10)$cluster
    			for (i in 1:n){ T[i,ind[i]] = 1 }
			TT[[i]] = T
 		}
 		V = switch(method,'svd'= fstep.fisher(Y,TT[[i]],kernel),'gs'= fstep.GramSc(Y,TT[[i]],kernel),'reg'  = fstep.qiao(Y,TT[[i]],kernel))
  		prms      = fem.mstep(Y,V,TT[[i]],model=model,method=method)
  		res.estep = fem.estep(prms,Y,V)
  		Ltmp[i]   = res.estep$loglik
	}
	T = TT[[which.max(Ltmp)]]
  }
   V = switch(method,'svd'= fstep.fisher(Y,T,kernel),'gs'= fstep.GramSc(Y,T,kernel),'reg'  = fstep.qiao(Y,T,kernel))
   prms = fem.mstep(Y,V,T,model=model,method=method)
   res.estep = fem.estep(prms,Y,V)
   Lobs[1] = res.estep$loglik
  
  # Main loop
  Linf_new  = Lobs[1]
  for (i in 1:maxit){
    # The three main steps F, M and E
    V = switch(method,
               'svd'= fstep.fisher(Y,T,kernel),
               'gs'= fstep.GramSc(Y,T,kernel),
               'reg'  = fstep.qiao(Y,T,kernel))
    prms      = fem.mstep(Y,V,T,model=model,method=method)
    if (prms$test !=0) stop("some classes become empty",call.=F)
    res.estep = fem.estep(prms,Y,V)
    T         = res.estep$T
    Lobs[i+1] = res.estep$loglik

    # Stop criterion
    if (i>=3){
      acc      = (Lobs[i+1] - Lobs[i]) / (Lobs[i] - Lobs[i-1])
      Linf_old = Linf_new
      Linf_new <- try( Lobs[i] + 1/(1-acc) * (Lobs[i+1] - Lobs[i]))
      if (is.na(Linf_new)) stop("some classes become empty\n",call.=F)
      if (abs(Linf_new - Linf_old) < eps) {break}
    } 
  }
  
  # Returning the results
  cls  = max.col(T)
  crit = fem.criteria(Lobs[(i+1)],T,prms,n)
  rownames(V) = colnames
  colnames(V) = paste('U',1:d,sep='')
  list(K=K,cls=cls,P=T,U=V,aic=crit$aic,mean=prms$mean,my=prms$my,prop=prms$prop,D=prms$D,model=prms$model,bic=crit$bic,icl=crit$icl,loglik=Lobs[2:(i+1)],ll=Lobs[i+1],method=method)
}

