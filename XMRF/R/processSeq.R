processSeq <-
function(X,quanNorm=.75,nLowCount=20,percentLowCount=.95,NumGenes=500,PercentGenes=.1)
{
  n = nrow(X)
  p = ncol(X)
  dis = apply(X,1,quantile,quanNorm)
  qnum = mean(dis)
  qX = X/(dis%o%rep(1,p)/qnum)    
  nums = apply(qX<=nLowCount,2,sum)
  indout = which(nums>(n*percentLowCount))
 
  # Wooi addition - 09/08/2014
  fX = qX
  if(length(indout) > 0){  fX = qX[,-indout] }
  
  tG = min(NumGenes,floor(p*PercentGenes))
  if(tG<ncol(fX)){
    lX = log(fX + 1)
    vars = apply(lX,2,var)
    or = order(vars,decreasing=TRUE)
    ffX = X[,or[1:tG]]
  }else{ffX = fX}  
  kslog = suppressWarnings(ks.test(c(log(ffX+1)),"ppois",mean(log(ffX+1)))$statistic)
  ks = NULL
  alphas = seq(.2,.5, l=30)      
  for(i in 1:length(alphas)){
    ks[i] = suppressWarnings(ks.test(c(ffX^alphas[i]),"ppois",mean(ffX^alphas[i]))$statistic)
  }
  if(kslog<min(ks)){
    Xout = log(ffX+1)
  }else{
    alphaopt = alphas[which.min(ks)]
    Xout = ffX^alphaopt
  }
 return(floor(Xout))
}
