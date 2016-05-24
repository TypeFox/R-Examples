llikeRPMMObject = function(o, x, type){
  J = dim(x)[2]
  L = rep(0, dim(x)[1])
  for(j in 1:J){
    if(type=="blc") {
      l = dbeta(x[,j], o$unsplit$a[j], o$unsplit$b[j], log=TRUE)
      L = L + ifelse(is.na(l),0,l)
    }
    else if(type=="glc"){
      l = dnorm(x[,j], o$unsplit$mu[j], o$unsplit$sigma[j], log=TRUE)
      L = L + ifelse(is.na(l),0,l)
    }
  }
  L
}

ebayes = function(rpmm,x,type,nodelist=NULL){
  ApplyFun = get(paste(type,"TreeApply",sep=""))
  if(is.null(nodelist)){
     nodelist = unlist(ApplyFun(rpmm,
      function(u,tree) u,
      terminalOnly=TRUE, asObject=FALSE))
  }
  NN = length(nodelist)
  nx = dim(x)[1]
  L = matrix(NA, nx, NN)
  eta = rep(NA, NN)
  for(i in 1:NN){
     L[,i] = llikeRPMMObject(rpmm[[nodelist[i]]],x,type)
     eta[i] = sum(rpmm[[nodelist[i]]]$weight)
  }
  eta = eta/sum(eta)

  L[is.na(L)] = min(L,na.rm=TRUE)
  for(i in 1:nx){
     L[i,] = exp(L[i,] - max(L[i,])) * eta
     L[i,] = L[i,]/sum(L[i,])
  }
  colnames(L) = nodelist
  L
}

glmLC = function(y,W,family=quasibinomial(),eps=1E-8,Z=NULL){
  K = dim(W)[2]
  C = as.vector(col(W))
  W = as.vector(W)
  y = rep(y,K)
 
  if(!is.null(Z)) {
    ZZ = NULL
    for(j in 1:dim(Z)[2]) {
      ZZ = cbind(ZZ,rep(Z[,j],K))
   }
   colnames(ZZ) = colnames(Z)
   ZZ=ZZ[W>eps,]
  }

  y=y[W>eps]
  C=C[W>eps]
  W=W[W>eps]
  C=as.factor(C)

  if(is.null(Z)) fit = glm(y~C, weights=as.vector(W), family=family)
  else fit = glm(y~C+ZZ, weights=as.vector(W), family=family)

  attr(fit,"Clevels") <- levels(C)
  fit
}

predict.glcTree = function(object, newdata=NULL, nodelist=NULL, type="weight", ...){
  x = object
  if(is.null(newdata)) {
    if(type=="weight") return(glcTreeLeafMatrix(x))
    else return(glcTreeLeafClasses(x))
  }

  w = ebayes(x, newdata, type="glc", nodelist=nodelist)  
  if(type=="weight") return(w)
  
  clnum = apply(w,1,which.max)
  as.factor(colnames(w))[clnum]
}

predict.blcTree = function(object, newdata=NULL, nodelist=NULL, type="weight", ...){
  x = object
  if(is.null(newdata)) {
    if(type=="weight") return(blcTreeLeafMatrix(x))
    else return(blcTreeLeafClasses(x))
  }

  w = ebayes(x, newdata, type="blc", nodelist=nodelist)  
  if(type=="weight") return(w)
  
  clnum = apply(w,1,which.max)
  as.factor(colnames(w))[clnum]
}

