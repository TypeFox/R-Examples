eigen.pda = function(pdaList,plotresult=TRUE,npts=501,...)
{
  rangval = pdaList$resfdlist[[1]]$basis$rangeval
  
  m = length(pdaList$resfdlist)  
  tfine = seq(rangval[1],rangval[2],length.out=npts)

  bwtlist = pdaList$bwtlist
  awtlist = pdaList$awtlist
  ufdlist = pdaList$ufdlist

  if(m == 1){
    d = length(bwtlist)
    xlabstr = names(bwtlist[[1]]$fd$fdnames)[[1]]

    betamat = array(0,c(npts,d,d))
      
    for(i in 1:d){
      betamat[,1,d-i+1] = -eval.fd(tfine,bwtlist[[i]]$fd)
      if(i < d) betamat[,i+1,i] = 1
    }
    
    if(!is.null(awtlist)){
      umat = matrix(0,npts,d)
      for(i in 1:length(awtlist)){
        umat[,1] = umat[,1] + 
          eval.fd(tfine,awtlist[[i]]$fd)*eval.fd(tfine,ufdlist[[i]])
      }
    
    }
    
  }
  else{
    d = length(bwtlist[[1]][[1]])
    xlabstr = names(bwtlist[[1]][[1]][[1]]$fd$fdnames)[[1]]
#    betamat = array(0,c(npts,m,m,d))

    betamat = array(0,c(npts,m*d,m*d))
    
    for(k in 1:d){
      for(j in 1:m){
        for(i in 1:m){
#                betamat[,i,j,k] = eval.fd(tfine,bwtlist[[i]][[j]][[k]]$fd)
          if(!is.null(bwtlist[[i]][[j]][[k]]))
          betamat[,j,m*(d-k)+i] = -eval.fd(tfine,bwtlist[[i]][[j]][[k]]$fd)
        }
        if(k < d){
            betamat[,m*k+j,m*(k-1)+j] = 1
        }
      }
    }
    
    if(!is.null(awtlist)){
      umat = matrix(0,npts,m*d)
      for(k in 1:d){
        for(i in 1:length(awtlist[[k]])){
          if(!is.null(awtlist[[k]][[i]]))
            umat[,k] = umat[,k] + 
               eval.fd(tfine,awtlist[[k]][[i]]$fd)*eval.fd(tfine,ufdlist[[k]][[i]])
        }
      }
    }
    
  } 
  
  eigvals = matrix(0,npts,m*d)
  limvals = matrix(0,npts,m*d)
              
  for(i in 1:npts){
      eigvals[i,] = eigen(betamat[i,,])$values
      if(!is.null(awtlist)) limvals[i,] = solve(betamat[i,,],umat[i,])
  }
  
  if(plotresult){
     if(!is.null(awtlist)) par(mfrow=c(3,1))
     else par(mfrow=c(2,1))
     matplot(tfine,Re(eigvals),type='l',xlab=xlabstr,ylab='Real',main='Eigenvalues',...)
     abline(h = 0)
     matplot(tfine,Im(eigvals),type='l',xlab=xlabstr,ylab='Imaginary',...)
     abline(h = 0)
    
     if(!is.null(awtlist))
     matplot(tfine,limvals[,1:d],type='l',xlab=xlabstr,ylab='Fixed Point',main='Instantaneous Limits',...)

  }
  
  return(list(argvals=tfine,eigvals=eigvals,limvals=limvals))  
} 
  
  