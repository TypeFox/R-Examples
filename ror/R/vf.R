sample.vfs.rejection <- function(performances, preferences,
                       nr=10000, updInterval=1000) {
  do.sample(performances, preferences, nr, thinning=1, updInterval, sampler=2)
}


sample.vfs.gibbs <- function(performances, preferences,
                       nr=10000, thinning=1, updInterval=1000) {
  do.sample(performances, preferences, nr, thinning, updInterval, sampler=1)  
}

do.sample <- function(performances, preferences, nr, thinning, updInterval, sampler) {    
  stopifnot(nr > 0)
  stopifnot(thinning > 0)
  
  ror <- vf.create(performances, nr, thinning, sampler)
  if (is.matrix(preferences)) {
    for (i in 1:nrow(preferences)) {
      ror.addPreference(ror, preferences[i,1], preferences[i,2])
    }
  }
  vf.sample(ror, updInterval)

  list(vfs=vf.allValueFunctions(ror, ncol(performances)),
       misses=vf.getMisses(ror))
}

vf.getMisses <- function(ror) {
  .jcall(ror$model, "I", method="getMisses")
}  

vf.allValueFunctions <- function(ror, nVf) {
  lapply(seq(1, nVf), function(x) {vf.getValueFunctionsForCriterion(ror, x)})
}

vf.create <- function(perfMat, nrVF, thinning, sampler) {
  model <- .jnew("fi/smaa/libror/r/ValueFunctionSamplerRFacade",
                 as.vector(perfMat), as.integer(nrow(perfMat)),
                 as.integer(nrVF), as.integer(thinning), as.integer(sampler))
  list(model=model, rownames=rownames(perfMat), colnames=colnames(perfMat))
}


vf.sample <- function(ror, upd) {
  .jcall(ror$model, "V", method="sample", as.integer(upd))
}

vf.getValueFunctionsForCriterion <- function(ror, cIndex) {
  .jcall(ror$model,
         "[[D",
         method="getValueFunctionsForCriterion",
         as.integer(cIndex-1),
         simplify=TRUE)
}


