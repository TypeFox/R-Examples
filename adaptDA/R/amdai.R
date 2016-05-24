amdai <-
function(X,cls,model='qda'){
  # Initialization
  C = max(cls)
  n = nrow(X)
  p = ncol(X)
  m = matrix(NA,C,p)
  prop = rep(c(NA),1,C)
  V = array(NA,c(C,p,p))
  
  # Learning
  for (i in 1:C){
    m[i,] = colMeans(X[cls==i,])
    prop[i] = nrow(X[cls==i,]) / n
    V[i,,] = cov(X[cls==i,])
  }
  
  if (model=='qda'){
    prms = list(model='qda',C=C,p=p,mean=m,prop=prop,var=V)
  }
  if (model=='lda'){
    VV = matrix(0,p,p)
    for (i in 1:C){
      VV = VV + prop[i] * V[i,,]
    }
    for (i in 1:C) {V[i,,] = VV}
    prms = list(model='lda',C=C,p=p,mean=m,prop=prop,var=V)
  }
  class(prms) <- "amdai"
  prms
}
