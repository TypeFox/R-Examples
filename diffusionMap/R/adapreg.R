adapreg = function(D,y,mmax=min(50,length(y)),fold=NULL,nfolds=10,nrep=5){

  # set up folds for CV
  if(length(fold)!=length(y)){ #defaults to nfolds-fold CV
    fold = sample(1:nfolds,length(y),replace=T)
  }    
  nfolds = length(table(fold))

  print(paste("Doing ",nfolds,"-fold cross-validation",sep=""))
  objtmp = Inf
  for(kk in 1:nrep){
    popt = optimize(adapreg.m,lower=0,upper=epsilonCompute(D,runif(1,.05,.2)),D=D,y=y,mmax=mmax,fold=fold, objfun=TRUE) # optimize CV risk over epsilon 

    if(popt$objective < objtmp){
      objtmp = popt$objective
      mintmp = popt$minimum
    }
  }
  ARopt = adapreg.m(mintmp,D,y,mmax=mmax,fold=fold)
  
  return(list(mincvrisk = objtmp, mopt = ARopt$mopt, epsopt = mintmp, y.hat= ARopt$y.hat,  coeff = ARopt$coeff))

}
