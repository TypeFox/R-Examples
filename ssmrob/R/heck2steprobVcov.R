heck2steprobVcov <-
function(y1vec, y2vec, x1Matr, x2Matr, eststage1, eststage2, eststage2sigma, weights = rep(1,nrow(y1vec)), t.c = 1.345)
{
  n=sum(y1vec==1)
  sumheck=0 #
  for(i in 1:n)  # sandwich estimator of the first term
  {
    sumheck=sumheck+PsiMest(x2Matr[i,],y2vec[i],eststage2,eststage2sigma,t.c,weights[i])%*%t(PsiMest(x2Matr[i,],y2vec[i],eststage2,eststage2sigma,t.c,weights[i]))
  }
  
  s1=matrix(0,length(eststage2)-1,length(eststage1$coeff)); 
  s2=0;   
  for(i in 1:n)
  {
    if((y2vec[i]-t(x2Matr[i,])%*%as.vector(eststage2))/eststage2sigma<(-t.c))
    {
      s1=s1+0
      s2=s2-t.c*weights[i]*drop(dLambdadSM(x1Matr[i,],eststage1$coeff))*(x1Matr[i,])
    }
    else if(((y2vec[i]-t(x2Matr[i,])%*%as.vector(eststage2)))/eststage2sigma>(t.c))
    {
      s1=s1+0
      s2=s2+t.c*weights[i]*drop(dLambdadSM(x1Matr[i,],eststage1$coeff))*(x1Matr[i,])
    }
    else {
      s1=s1+( x2Matr[i,1:dim(x2Matr)[2]-1] )%*%t(x1Matr[i,])*eststage2[dim(x2Matr)[2]]*weights[i]*drop(dLambdadSM(x1Matr[i,],eststage1$coeff))/eststage2sigma
      s2=s2+(x2Matr[i,dim(x2Matr)[2]]*eststage2[dim(x2Matr)[2]])*weights[i]*dLambdadSM(x1Matr[i,],eststage1$coeff)/eststage2sigma*x1Matr[i,]
    }
  }
  xdx=rbind(s1,s2)
  term2=xdx%*%vcov(eststage1)%*%t(xdx)  
  
  result=(solve(MmatrM(x2Matr,y2vec,eststage2,eststage2sigma,t.c,weights))%*%   # asymptotic variance for Heckman-M-estimator
    (sumheck+term2)%*%solve(MmatrM(x2Matr,y2vec,eststage2,eststage2sigma,t.c,weights)))
  return(result)
}
