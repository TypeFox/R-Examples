soft <- function(X, knowns, P=NULL, k=ifelse(!is.null(P),ncol(P),ifelse(!is.null(B),ncol(B),length(unique(class)))), B=NULL, class=NULL, init.params=init.model.params(X, knowns, class=class, B=P, k=k), model.structure=getModelStructure(), stop.likelihood.change=10^-5, stop.max.nsteps=100, trace=FALSE, b.min=0.025, all.possible.permutations=FALSE) {
  if (is.null(dim(knowns)) || is.data.frame(knowns)) knowns = as.matrix(knowns)
  if (is.null(dim(X)) || is.data.frame(X)) X = as.matrix(X)
  if (is.null(P)) {
    if (!is.null(B)) {
      P=B
      if (is.null(k) | k < 2)
          k=ncol(P)
    } else if (!is.null(class)){
      P = get.simple.beliefs(class, b.min=b.min)
      if (is.null(k) | k < 2)
          k = length(unique(class))
    } else 
      stop("Argument P need to be specified")
  }
  if (k > ncol(P))
    P = cbind(P, matrix(0,nrow(P),k - ncol(P)))
  if (ncol(X) != ncol(knowns))  
      stop("number of columns in X and knowns must agree")
  init.params$P = P
  init.params$m = nrow(knowns)
  init.params$n = nrow(knowns) + nrow(X)
  init.params$k = k
  init.params$d = ncol(X)
  result = bgmm.internal(internal.funct=soft.internal, X=rbind(knowns, X), init.params=init.params, model.structure=model.structure, stop.likelihood.change=stop.likelihood.change, stop.max.nsteps=stop.max.nsteps, trace=trace, all.possible.permutations=all.possible.permutations)

  result$likelihood = result$likelihood           #loglikelihood.mModel(result, X)
  result$X = X
  result$knowns = knowns
  result$model.structure = model.structure
  result$B = B
  class(result) = c("softModel", "mModel")
  if (!is.null(colnames(X))) {
      dimnames(result$cvar) = list(NULL, colnames(X), colnames(X))
  }
  if (!is.null(colnames(knowns))) {
      dimnames(result$cvar) = list(NULL, colnames(knowns), colnames(knowns))
  }
  
  result$dof = getDFinternal(result)

  result
}


# bgmm.internal.call

soft.internal <- function(X, model.params, model.structure, stop.likelihood.change=10^-5, stop.max.nsteps=100, trace=F) {
  prev.likelihood = -Inf
  n.steps = 0
  repeat {
    n.steps = n.steps +1
    tmp = soft.e.step(X, model.params) 
    model.params = bgmm.m.step(X, model.params, model.structure, tmp$tij, priors.like.bgmm=FALSE)
    if (trace) {
      cat("step:          ", n.steps, "\n likelihood:   ", tmp$log.likelihood, "\n change:       ", tmp$log.likelihood - prev.likelihood, "\n\n")
    }
    if ((abs(tmp$log.likelihood - prev.likelihood)/ifelse(is.infinite(prev.likelihood), 1,  (1+abs(prev.likelihood))) < stop.likelihood.change) || 
        (n.steps >= stop.max.nsteps)) {
          model.params$likelihood = tmp$log.likelihood
          break
      }
    prev.likelihood = tmp$log.likelihood
  }
  model.params$likelihood = prev.likelihood
  model.params$n.steps = n.steps
  model.params$tij = tmp$tij
  model.params
}

unsupervised <- function(X, k, init.params=init.model.params(X, knowns=NULL, k=k), model.structure=getModelStructure(), stop.likelihood.change=10^-5, stop.max.nsteps=100, trace=FALSE) {
  if (is.null(dim(X)) || is.data.frame(X)) X = as.matrix(X)
  init.params$P = NULL
  init.params$m = 0
  init.params$n = nrow(X)
  init.params$k = k
  init.params$d = ncol(X)
  result = soft.internal(X, init.params, model.structure, stop.likelihood.change=stop.likelihood.change, stop.max.nsteps=stop.max.nsteps, trace=trace)
  result$likelihood = result$likelihood           #loglikelihood.mModel(result, X)
  result$X = X
  result$knowns = NULL
  result$model.structure = model.structure
  result$B = NULL
  class(result) = c("softModel", "mModel")
  if (!is.null(colnames(X))) {
      dimnames(result$cvar) = list(NULL, colnames(X), colnames(X))
  }

  result$dof = getDFinternal(result)

  result
}

