mModelList <- function(X, knowns, B=NULL, P=NULL, class=NULL, kList=ncol(B), init.params=NULL, stop.likelihood.change=10^-5, stop.max.nsteps=100, trace=FALSE, mean = c("D","E"), between = c("D","E"), within = c("D","E"), cov = c("D","0"), funct=belief, all.possible.permutations=FALSE) {
  modelList = NULL
  nazw = NULL
  loglike = NULL
  params = NULL
  ind = 0
  d = ncol(X)
  # 
  # for d=1 not all combinations have meaning
  if (d<2) {
    within = "D"
    cov = "D"
  }
  
  for (k in kList) {
    if (is.null(init.params)) {
        model.paramsK = init.model.params(X, knowns, B=B, P=P, class=class, k=k)
      } else {
        model.paramsK = init.params
      }

      for (i in mean) {
        for (j in between) {
          for (m in within) {
            for (l in cov) {
              model.settings = list(mean = i, between = j, within = m, cov = l)
              if (trace) 
                 cat("--------------------------------\n",nazw[[ind]],"\n")
                ind = ind + 1
                nazw[[ind]] = paste("k=",k,"  structure=",i,j,m,l,sep="")
                
                modelList[[ind]] = funct(X, knowns, B=B, P=P, class=class, k=k, init.params=model.paramsK, model.structure=model.settings, stop.max.nsteps=stop.max.nsteps, trace=trace, all.possible.permutations=all.possible.permutations)
                loglike[[ind]] = modelList[[ind]]$likelihood 
                params[[ind]] = ifelse(i=="D", k*d, d) + # mean value
                                  ifelse(j=="D", k, 1) * # variance between clusters
                                  (ifelse(m=="D", d, 1)+ # diagonal
                                   ifelse(l=="0", 0, 1)*ifelse(m=="D", d*(d-1)/2, 1) ) # covariances
              }
            }
          }
    }
  }
  result = list(models=modelList, loglikelihoods=loglike, names=nazw, params=params, kList = kList)
  class(result) = "mModelList"
  result
}


beliefList <- function(..., funct=belief)  mModelList(..., funct=belief)

semisupervisedList <- function(..., funct=semisupervised)  mModelList(..., funct=semisupervised)

softList <- function(..., funct=soft)  mModelList(..., funct=soft)

unsupervisedList <- function(X, kList = 2, ...) {
   mModelList(X=X[-1,,drop=F], knowns=X[1,,drop=F], B=matrix(1/min(kList),1,min(kList)), ..., funct=soft)
}

