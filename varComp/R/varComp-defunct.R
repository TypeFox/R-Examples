varComp.LinScore.LC12 <-
function(null.fit, X, Y, K, null, w, ...)
{ 
#Actual testing function for LC12. 
#i.	null.fit: object from varComp under the null. 
#ii.	X: the same as in varComp.test.
#iii.	Y: the same as in varComp.test. 
#iv.	K: the same as in varComp.test, after pre-/post-multiplying Z matrices, if applicable. 
#v.	null: the same as in varComp.test. 
#vi.	w: weights of scores (not used for the LC12 method). 
  if(!all(w==1)) warning('Weights are ignored in LC12 approximation method')  ## FIXME
  n=nrow(K[[1]])
  nNull=length(null)
  nK=length(K)
  if(nK!=1L+nNull){ # FIXME: implement
    return(structure(list(p.value=NA_real_),class='htest'))
  }
  sigma=diag(null.fit$sigma2,n)+Reduce('+', mapply('*', null.fit$varComps, K[null], SIMPLIFY=FALSE))
  sigma.ihalf=local({eg=eigen(sigma,TRUE); tcrossprod( sweep(eg$vec, 2L, sqrt(eg$val), '/'), eg$vec)})
  
  qrix=qr(sigma.ihalf%*%X)
  P01=sigma.ihalf%*%(diag(1,n)-tcrossprod(qr.Q(qrix)[,seq_len(qrix$rank),drop=FALSE]))%*%sigma.ihalf
  
  K3=K[-null][[1]]  ## FIXME: Check K3=Reduce('+', K[-null])
  SI=.5*crossprod(Y, P01%*%K3%*%P01%*%Y)

  Delta=matrix(NA_real_, nK, nK)
  Phi=numeric(nK)
    null1=c(null, nK+1)
    K[[nK+1]]=diag(1, n)
  for(i in 1:nK) {
    for(j in i:nK) Delta[i,j]=Delta[j,i]=sum(diag(P01%*%K[[null1[i]]]%*%P01%*%K[[null1[j]]]))
    Phi[i]=sum(diag(P01%*%K3%*%P01%*%K[[null1[i]]]))
  }
  nuI=.5*sum(diag(P01%*%K3%*%P01%*%K3))-.5*crossprod(Phi, solve(Delta, Phi))
  deltaI=.5*sum(diag(P01%*%K3))
  aI=nuI/2/deltaI
  gI=2*deltaI*deltaI/nuI
  pval=pchisq(SI/aI, gI, lower.tail=FALSE)
  ans=list(statistic=c(SI=SI), 
           p.value=pval, 
           alternative='greater', 
           parameter=c(scale=aI, df=gI), 
           null.value=structure(0, names='interaction variance component'), null.fit=null.fit, 
           method='Li and Cui (2012) Satterwaite-approximated Variance Component Test'
          )
  class(ans)='htest'
  ans
}

varComp.LinScore.LC12Boundary <-
function(null.fit, X, Y, K, null, ...)
{
#Boundary corrected LC12 method. Argument definitions are the same as in varComp.LinScore.LC12. 
  .Deprecated('varComp.test.aproximation.Satterthwaite')
  n=nrow(K[[1]])
  nNull=length(null)
  nK=length(K)
  stopifnot(nK==1L+nNull) ## FIXME
  sigma=diag(null.fit$sigma2,n)+Reduce('+', mapply('*', null.fit$varComps, K[null], SIMPLIFY=FALSE))
  sigma.ihalf=local({eg=eigen(sigma,TRUE); tcrossprod( sweep(eg$vec, 2L, sqrt(eg$val), '/'), eg$vec)})
  
  qrix=qr(sigma.ihalf%*%X)
  P01=sigma.ihalf%*%(diag(1,n)-tcrossprod(qr.Q(qrix)[,seq_len(qrix$rank),drop=FALSE]))%*%sigma.ihalf
  
  K3=K[-null][[1]]  ## FIXME: Check K3=Reduce('+', K[-null])
  SI=.5*crossprod(Y, P01%*%K3%*%P01%*%Y)

  Delta=matrix(NA_real_, nK, nK)
  Phi=mu12=S12=numeric(nK)
    null1=c(null, nK+1)
    K[[nK+1]]=diag(1, n)
  for(i in 1:nK) {
    for(j in i:nK) Delta[i,j]=Delta[j,i]=sum(diag(P01%*%K[[null1[i]]]%*%P01%*%K[[null1[j]]]))
    Phi[i]=sum(diag(P01%*%K3%*%P01%*%K[[null1[i]]]))
    mu12[i]=.5*sum(diag(P01%*%K[[null1[i]]]))
    S12[i]=.5*crossprod(Y, P01%*%K[[null1[i]]]%*%P01%*%Y)
  }
  nuI=.5*sum(diag(P01%*%K3%*%P01%*%K3))-.5*crossprod(Phi, solve(Delta, Phi))
  deltaI=.5*sum(diag(P01%*%K3))-crossprod(Phi, solve(Delta, mu12-S12))
  aI=nuI/2/deltaI
  gI=2*deltaI*deltaI/nuI
  pval=pchisq(SI/aI, gI, lower.tail=FALSE)
  ans=list(statistic=c(SI=SI), 
           p.value=pval, 
           alternative='greater', 
           parameter=c(scale=aI, df=gI), 
           null.value=structure(0, names='interaction variance component'), null.fit=null.fit, 
           method='Corrected Li and Cui (2012) Satterwaite-approximated Variance Component Test'
          )
  class(ans)='htest'
  ans
}


varComp.LinScore.approximation <-
function(approximation, ...)
{
#A Switcher that choose the approximation method for LinScore test when the null has additional variance components other than the error variance. 
#i.	approximation: A character scalar, specifying the mthod of approximation. 
#ii.	? Other arguments to be passed to the actual approximate testing function. 
  this.call=match.call()
  idx=which(names(this.call)=='approximation')
  stopifnot(length(idx)>0L && length(approximation)==1L)
  this.call=this.call[-idx]
  #this.call[[1]]=as.name(paste('varComp.LinScore', approximation, sep='.'))

  this.env=sys.frame(sys.parent(0L))
  parent.env(this.env)=parent.frame()

  fname=paste('varComp.LinScore', approximation, sep='.')
  toCall=get(fname, mode='function')
  eval(substitute({environment(zzz)=this.env}, list(zzz=as.name(fname))))

  this.call[[1]]=as.name(fname)
  eval(this.call)
  
#  idx=which(names(this.call)=='envir')
#  if(length(idx)>0L)this.call=this.call[-idx]
#  eval(this.call, envir=parent.frame()) # explicit envir=parent.frame() is needed!
#  this.env=sys.frame(sys.parent(0L))
#  bak.parent=parent.env(this.env)
#  parent.env(this.env)=parent.frame()
#  ans=eval(this.call)
#  parent.env(this.env)=bak.parent
#  ans
}

