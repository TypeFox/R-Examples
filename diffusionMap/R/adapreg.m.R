adapreg.m = function(epsilon,D,y,mmax=min(50,length(y)),fold=NULL,nfolds=10,objfun=FALSE){

  print(paste("Finding optimal model for epsilon=",round(epsilon,4)))
  
  dmap = diffuse(as.matrix(D),eps.val=epsilon,neigen=mmax,t=1)
  
  pred.y = matrix(0,length(y),mmax)

  if(length(fold)!=length(y)){ #defaults to nfolds-fold CV
    fold = sample(1:nfolds,length(y),replace=T)
  }    
  nfolds = length(table(fold))
    
  for(ii in 1:nfolds){
    if(sum(fold==ii)>1){
        # regress y on diffusion map coordinates
      AR = lm(y[fold!=ii]~dmap$X[fold!=ii,])
      for(jj in 1:mmax){ # predictions for each size model
        pred.y[fold==ii,jj] = cbind(rep(1,sum(fold==ii)),dmap$X[fold==ii,1:jj])%*%AR$coeff[1:(jj+1)]
      }
    }
    if(sum(fold==ii)==1){ # if only 1 data point in the fold
      # regress y on diffusion map coordinates
      AR = lm(y[fold!=ii]~dmap$X[fold!=ii,])
      for(jj in 1:mmax){ # predictions for each size model
        pred.y[fold==ii,jj] = c(1,dmap$X[fold==ii,1:jj])%*%AR$coeff[1:(jj+1)]                
      }
    }
  }

  # compute MSE
  risk.dmap = apply((pred.y-y)^2,2,mean)
  mopt = which.min(risk.dmap)

  ARopt = lm(y~dmap$X[,1:mopt])
        
  if(objfun==TRUE){
    return(min(risk.dmap))
  } else{
    return(list(mincvrisk = min(risk.dmap), mopt = mopt, cvrisks = risk.dmap, epsilon = epsilon, y.hat = ARopt$fitted.values, coeff = ARopt$coeff))
  }
}
