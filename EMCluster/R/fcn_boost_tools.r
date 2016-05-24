### This file contains tool functions for boosting.
### Written: Wei-Chen Chen on 2008/11/02.


### This function returns an initial result in a list
###   model: list[nclasses]
###   result: data.frame[nclasses, 2]
em.init <- function(x, nclasses, method.init = c("emgroup", "em.EM")){
  if(method.init[1] == "emgroup"){
#    model <- lapply(nclasses, function(i){ emgroup(x, nclass = i) })
    x <- t(x)
    model <- lapply(nclasses, function(i){ emgroup.wt(x, nclass = i) })
    model <- lapply(model, wt2wot)
  } else if(method.init[1] == "em.EM"){
    model <- lapply(nclasses, function(i){ em.EM(x, nclass = i) })
    ### WCC: TBD.
    # x <- t(x)
    # model <- lapply(nclasses, function(i){ em.EM.wt(x, nclass = i) })
    # model <- lapply(model, wt2wot)
  } else{
    stop("The initial method is not found.")
  }

  result <- data.frame(nclass = nclasses, bic = sapply(model, em.bic))

  ret <- list(model = model, result = result)
  class(ret) <- "initret"
  ret
}

print.initret <- function(x, ...){
  print(x$result)
}


em.Model <- function(x, initobj){
  nclass <- which.min(initobj$result$bic)
  ret <- emcluster(x, initobj$model[[nclass]])
  ret
}

