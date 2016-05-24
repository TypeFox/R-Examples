predict.amdai <-
function(object,Y,K,maxit=25,eps=1e-5,...){
  # Usage: res <- amda(prms,Y,K)
  # Initialization
  C = object$C
  p = object$p
  m = object$mean
  prop = object$prop
  V = object$var
  model = object$model
  ny = nrow(Y)
  
  # New objects
  L = rep(c(-Inf),1,(maxit+1))
  
  # Case K < C 
  if (K < C) stop('amda: The number of components must be higher or equal to ',C,'!')
  
  # Case K = C 
  if (K==C){
    #cat('amda: Usual classification (K=C)\n')
    T = .amdai.estep(object,Y,C,K)
    cls = max.col(T)
    prms = object
    prms$K = K
    i = 1;
    L[i+1] = .amdai.loglikelihood(prms,T,Y,C,K)
  }
  
  # Case K > C 
  if (K > C){
    #cat('amda: Classification with unobserved classes (K>C)\n')
    prms = .amdai.initparam(object,Y,K)
    for (i in 1:maxit){
      #cat('.')
      T = .amdai.estep(prms,Y,C,K)
      prms = .amdai.mstep(object,T,Y,C,K)
      if (sum(prms$prop < 1e-3) == 1){
        cat('!!! Empty class !!!\n')
        L[i+1] = -Inf
        break
      }
      L[i+1] = .amdai.loglikelihood(prms,T,Y,C,K)
      if (abs(L[i+1]-L[i])<ny*eps) break
    }
    #cat('\n')
  }
  # Returning the results
  cls = max.col(T)
  crit = .amdai.crit(L[(i+1)],T,prms,ny);
  res = list(model=model,C=prms$K,p=p,mean=prms$mean,prop=prms$prop,var=prms$var,cls=cls,P=T,crit=crit,loglik=L[2:(i+1)])
  class(res) <- "amdai"
  res
}
