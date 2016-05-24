predict.mModel <- function(object, X, knowns=NULL, B=NULL, P=NULL, ...) {
  # densities for all components
  if (is.null(dim(X)) || is.data.frame(X)) X = as.matrix(X)
  if (!is.null(knowns)) {
      if (is.null(dim(knowns)) || is.data.frame(knowns)) knowns = as.matrix(knowns)
      X = rbind(knowns, X)         # change
  }
  if ((is.null(B) & is.null(P) & !is.null(knowns)) | ({!is.null(B) | !is.null(P)} & is.null(knowns))) {
     stop("If knowns are specified there should be also B or P specified as well!")
  }
  lfik <- matrix(0, nrow(X), object$k)
  rownames(lfik) = rownames(X)
  for (i in 1:object$k) {
      if (object$d > 1) {
         ss = svd(object$cvar[i,,])
         rtas <- ss$d
         matc = t(ss$u[rtas > 10^-8, ]) %*% diag(rtas[rtas > 10^-8]^(-1/2)) %*% ss$v[rtas > 10^-8,]
         tx = apply(X, 1, get("-"), object$mu[i,,drop=F])
#         lfik[,i] <- exp(-colSums((matc %*% tx)^2)/2)/sqrt(prod(2*pi*rtas[rtas > 10^-8]))
         lfik[,i] <- -colSums((matc %*% tx)^2)/2 -sum(log(2*pi*rtas[rtas > 10^-8]))/2
      } else {
         tx = apply(X, 1, get("-"), object$mu[i,,drop=F])
         lfik[,i] <- -(tx^2)/(2*object$cvar[i,,]) - log(2*pi*object$cvar[i,,])/2
       }
  }
  
  fik <- exp(lfik)
  # numeric problems, trying to adjust, by keeping likelihood not so far from each other
  if (min(apply(fik, 1, max)) == 0) {
    lfik <- t(apply(lfik, 1, function(x) x - max(x)))
    fik = exp(lfik)
  }

  b.pi <- repeat.rows(object$pi, nrow(X))
  if (!is.null(B))          # change
        b.pi[1:nrow(knowns),] = B
#        b.pi[nrow(X) - (nrow(knowns):1) +1,] = B
  if (!is.null(P))          # change 
        b.pi[1:nrow(knowns),] = P * b.pi[1:nrow(knowns),]
#        b.pi[nrow(X) - (nrow(knowns):1) +1,] = P * b.pi[nrow(X) - (nrow(knowns):0),]

  tij =  t(apply(fik * b.pi, 1, normalize))
  class = get.labels.from.beliefs(tij)
  # return predictions as two separate slots
  tij.knowns = NULL
  tij.X = tij
  class.knowns = NULL
  class.X = class
  if (!is.null(P) | !is.null(B)) {
     tij.knowns = tij[1:nrow(knowns),,drop=F]
     tij.X = tij[-(1:nrow(knowns)),,drop=F]
     class.knowns = class[1:nrow(knowns)]
     class.X = class[-(1:nrow(knowns))]
  }
  list(tij.X=tij.X, tij.knowns = tij.knowns, class.X=class.X, class.knowns=class.knowns)
}


supervised <- function(knowns, class=NULL, k=length(unique(class)), B=NULL, P=NULL, model.structure=getModelStructure(), ...) {
  if (is.null(dim(knowns)) || is.data.frame(knowns)) knowns = as.matrix(knowns)
  if (is.null(class)) {
    if (!is.null(B)) {
      class = apply(B, 1, which.max)
      if (is.null(k) | k < 2)
          k=length(unique(class))
    } else if (!is.null(P)){
      class = apply(P, 1, which.max)
      if (is.null(k) | k < 2)
          k=length(unique(class))
    } else 
      stop("Argument class need to be specified")
  }
  result = init.model.params.knowns(knowns, class, k, ncol(knowns)) 

  # new means 
  if  (model.structure$mean!="D") {
    result$mu = repeat.rows(colMeans(knowns), k)
  }
  # are variances equal?
  if (model.structure$between=="E") {
      # averaging among clusters
     ncvar = matrix(0, result$d, result$d)
     for (i in 1:k) 
        ncvar = ncvar + result$cvar[i, , ] * result$pi[i]
     for (i in 1:k) 
        result$cvar[i, , ] = ncvar
  }
  if (model.structure$within=="E" && ncol(knowns)>1) {
      # averaging among variables
     for (i in 1:k) {
        ndiag = sum(diag(result$cvar[i, , ]))
        sdiag = ndiag/(ncol(knowns))
        result$cvar[i, , ] = min(sdiag, (sum(result$cvar[i, , ])-ndiag)/(ncol(knowns)*(ncol(knowns)-1)))
        diag(result$cvar[i, , ]) = sdiag
     }
  }
  # are covariances equal to 0?
  if (model.structure$cov=="0") {
   # covariances are equal to 0
   for (i in 1:k) 
      result$cvar[i, , ] = diag(diag(result$cvar[i, , ]), nrow=ncol(knowns))
  }  
  result$m = nrow(knowns)
  result$n = nrow(knowns)
  result$k = k
  result$d = ncol(knowns)
  result$knowns = knowns
  result$class = class
  class(result) = c("supervisedModel", "mModel")
  if (!is.null(colnames(knowns))) {
      dimnames(result$cvar) = list(NULL, colnames(knowns), colnames(knowns))
  }
  
  result$dof = getDFinternal(result)
  
  result
}

semisupervised <- function(X, knowns, class=NULL, k=ifelse(!is.null(class),length(unique(class)),ifelse(!is.null(B),ncol(B),ncol(P))),B=NULL,P=NULL, ...,  all.possible.permutations=FALSE) {
  if (is.null(dim(knowns)) || is.data.frame(knowns)) knowns = as.matrix(knowns)
  if (is.null(dim(X)) || is.data.frame(X)) X = as.matrix(X)
  if (is.null(class)) {
    if (!is.null(B)) {
      class = apply(B, 1, which.max)
      if (is.null(k) | k < 2)
          k=length(unique(class))
    } else if (!is.null(P)){
      class = apply(B, 1, which.max)
      if (is.null(k) | k < 2)
          k=length(unique(class))
    } else 
      stop("Argument class need to be specified")
  }
  if (ncol(X) != ncol(knowns))  
      stop("number of columns in X and knowns must agree")
  result = soft(X, knowns, get.simple.beliefs(class, k=k, b.min=0), k=k, ..., all.possible.permutations=all.possible.permutations) 
  result$X = X
  result$knowns = knowns
  result$class = class
  class(result) = c("semisupervisedModel", "mModel")
  result
}

belief <- function(X, knowns, B=NULL, k=ifelse(!is.null(B),ncol(B),ifelse(!is.null(P),ncol(P),length(unique(class)))), P=NULL, class=map(B), init.params=init.model.params(X, knowns, B=B, P=P, class=class, k=k), model.structure=getModelStructure(), stop.likelihood.change=10^-5, stop.max.nsteps=100, trace=FALSE, b.min=0.025,  all.possible.permutations=FALSE) {
  if (is.null(dim(knowns)) || is.data.frame(knowns)) knowns = as.matrix(knowns)
  if (is.null(dim(X)) || is.data.frame(X)) X = as.matrix(X)
  if (is.null(B)) {
    if (!is.null(P)) {
      B=P
      if (is.null(k) | k < 2)
          k=ncol(B)
    } else if (!is.null(class)){
      B = get.simple.beliefs(class, b.min=b.min)
      if (is.null(k) | k < 2)
          k=length(unique(class))
    } else 
      stop("Argument B need to be specified")
  }
  if (k > ncol(B))
    B = cbind(B, matrix(0,nrow(B),k - ncol(B)))
  if (ncol(X) != ncol(knowns))  
      stop("number of columns in X and knowns must agree")
  init.params$B = B
  init.params$m = nrow(knowns)
  init.params$n = nrow(knowns) + nrow(X)
  init.params$k = k
  init.params$d = ncol(X)
  result = bgmm.internal(internal.funct=belief.internal, X=rbind(knowns, X), init.params=init.params, model.structure=model.structure, stop.likelihood.change=stop.likelihood.change, stop.max.nsteps=stop.max.nsteps, trace=trace, all.possible.permutations=all.possible.permutations)
  result$X = X
  result$knowns = knowns
  result$B = B
  result$model.structure = model.structure
  if (!is.null(colnames(X))) {
      dimnames(result$cvar) = list(NULL, colnames(X), colnames(X))
  }
  if (!is.null(colnames(knowns))) {
      dimnames(result$cvar) = list(NULL, colnames(knowns), colnames(knowns))
  }

  result$dof = getDFinternal(result)

  class(result) = c("beliefModel", "mModel")
  result
}

#
# do we need to consider all possible permutations?
#
bgmm.internal <- function(internal.funct=belief.internal, init.params, ..., all.possible.permutations=all.possible.permutations) {
  if (all.possible.permutations) {
    resmax = NULL
    likemax = -Inf
    lpermut = permn(seq_along(init.params$mu))
    for (perm in lpermut) {
       tmodel.params = init.params
       tmodel.params$mu = tmodel.params$mu[perm,,drop=F]
       tmodel.params$cvar = tmodel.params$cvar[perm,,,drop=F]
       tmodel.params$pi = tmodel.params$pi[perm,drop=F]
       tmpr = internal.funct(..., model.params=tmodel.params)
       if (tmpr$likelihood > likemax) {
            likemax = tmpr$likelihood
            resmax = tmpr
       }
    }
    return(resmax)
  } 
  return(internal.funct(..., model.params=init.params))
}
 

# bgmm.internal.call

belief.internal <- function(X, model.params, model.structure, stop.likelihood.change=10^-5, stop.max.nsteps=100, trace=F) {
  prev.likelihood = -Inf
  n.steps = 0
  stopP = FALSE
  repeat {
    n.steps = n.steps +1
    tmp = bgmm.e.step(X, model.params) 
    model.params = bgmm.m.step(X, model.params, model.structure, tmp$tij)
    if (stopP)
          break
    if (trace) {
      cat("step:          ", n.steps, "\n likelihood:   ", tmp$log.likelihood, "\n change:       ", tmp$log.likelihood - prev.likelihood, "\n\n")
    }
    if ((abs(tmp$log.likelihood - prev.likelihood)/ifelse(is.infinite(prev.likelihood), 1,  (1+abs(prev.likelihood))) < stop.likelihood.change) || 
        (n.steps >= stop.max.nsteps)) {
          model.params$likelihood = tmp$log.likelihood
          stopP = TRUE
      }
    prev.likelihood = tmp$log.likelihood
  }
  
  model.params$likelihood = prev.likelihood
  model.params$n.steps = n.steps
  model.params$tij = tmp$tij
  model.params
}

#
# usefull functions
#

normalize <- function(x) x/sum(x)

repeat.rows <- function(x, k) matrix(x, k, length(x), byrow=T)

get.labels.from.beliefs <- function(B) {apply(B, 1, function(x) which.max(x)[1])}
map = get.labels.from.beliefs


determinant.numeric <- function (x, logarithm = TRUE, ...) {list(modulus = ifelse(logarithm,log(x),x), sign=sign(x))}

fij.mModel <- function(model, X) {
  # densities for all components
  lfik <- matrix(0, nrow(X), model$k)
  for (i in 1:model$k) {
      if (model$d > 1) {
         ss = svd(model$cvar[i,,])
         rtas <- ss$d
         matc = t(ss$u[rtas > 10^-8, ]) %*% diag(rtas[rtas > 10^-8]^(-1/2)) %*% ss$v[rtas > 10^-8,]
         tx = apply(X, 1, get("-"), model$mu[i,,drop=F])
         lfik[,i] <-  -colSums((matc %*% tx)^2)/2 - sum(log(2*pi*rtas[rtas > 10^-8]))/2
       } else {
        lfik[,i] <-   dnorm(X, model$mu[i,,drop=F], sqrt(model$cvar[i,,]), log=T)
       }
  }
  exp(lfik + repeat.rows(log(model$pi), nrow(X)))
#  fb.ik
}

loglikelihood.mModel <- function(model, X) {
  sum(log(rowSums(fij.mModel(model, X))))
}


#
# ICs
#

getGIC <- function(model, p=2, whichobs="unlabeled") {
  tmpDF = getDF(model)
  workingX = switch (whichobs, 
              unlabeled = model$X,
              labeled   = model$knowns,
              all       = rbind(model$X, model$knowns),
              stop("getGIC: unsupported value of the parameter whichobs"))
  n = max(nrow(workingX),1)

  fij   = fij.mModel(model, workingX)
  if ("numeric" %in% class(p)) return (-2*sum(log(rowSums(fij))) +  p * tmpDF)
  penalty = 0
  if ("character" %in% class(p)) 
      penalty = switch(p, 
          AIC     =  2 * tmpDF,
          AIC3    =  3 * tmpDF,
          AIC4    =  4 * tmpDF,
          AICc    = 2 * tmpDF + 2 * tmpDF * (tmpDF + 1) / (max(n - tmpDF - 1,1)),
          AICu    = 2 * tmpDF + 2 * tmpDF * (tmpDF + 1) / (max(n - tmpDF - 1,1)) + n * log(n/max(n - tmpDF - 1,1)),
          CAIC    = tmpDF * (1 + log(n)),
          BIC     = log(n) * tmpDF,
          MDL     = log(n) * tmpDF,
          CLC     = 2 * getEntropy(fij), 
          ICLBIC  = log(n) * tmpDF + 2 * getEntropy(fij),
          AWE     = 2 * tmpDF * (3/2 + log(n)),
          stop("getGIC: unkown penalty description"))
  
 -2*sum(log(rowSums(fij))) + penalty
}

# degress of fredom
getDFinternal <- function(model) {
  k = model$k
  d = model$d
  ifelse(model$model.structure$mean=="D", k*d, d) + # mean value
                    ifelse(model$model.structure$between=="D", k, 1) * # variance between clusters
                    (ifelse(model$model.structure$within=="D", d, 1)+ # diagonal
                     ifelse(model$model.structure$cov=="0", 0, 1)*ifelse(model$model.structure$between=="D", d*(d-1)/2, 1) ) # covariances
}

getDF <- function(model) {
  if (is.null(model$dof)) {
      return(getDFinternal(model))
  } 
  model$dof
}

getEntropy <- function(matr) {
  matr <- matr*log(matr)
  matr[is.nan(matr)] = 0
  -sum(matr)
}

chooseModels <- function(models, kList = NULL, struct = NULL) {
   if (!is.null(struct)) {  
    indStr = unlist(lapply(struct, function(x) {grep(x, models$names)}))
   } else {
    indStr = seq_along(models$names)
   }
   if (!is.null(kList)) {  
    indList = unlist(lapply(paste("k=",kList," ",sep=""), function(x) {grep(x, models$names)}))
   } else {
    indList = seq_along(models$names)
    kList = models$kList
   }
   indLS = intersect(indList, indStr)
   models2 = models
   models2$kList = kList
   models2$loglikelihoods = models$loglikelihoods[indLS]
   models2$names = models$names[indLS]
   models2$params = models$params[indLS]
   models2$models = models$models[indLS]
   models2
}



chooseOptimal <- function(models, penalty=2, ...) {
   values = sapply(models$models, getGIC, p=penalty, ...)
   models$models[[which.min(values)[1]]] 
}



crossval <- function(model=NULL, X=NULL, knowns=NULL, class=NULL, k=length(unique(class)),B=NULL,P=NULL, model.structure=getModelStructure(), ..., folds = 2, fun=belief) {
   if (!is.null(model)) {
      X = model$X
      knowns = model$knowns
      class = model$class
      B = model$B
      P = model$P
      k = model$k
      model.structure = model$model.structure
   }
   if (is.null(knowns)) {
      stop("Argument knowns has to be specified")
   }
   if (!is.data.frame(knowns) && !is.matrix(knowns)) {
      stop("Argument knowns has to be of the class matrix or data frame")
   }
   if (nrow(knowns) < folds) {
      stop("Number of folds cannot be greater than number of known cases")
   }
   m = nrow(knowns)
   n = nrow(X)
   indKnowns = sample(rep(1:folds, length.out=m))
   indX      = sample(rep(1:folds, length.out=n))
   errors = numeric(folds)
   for (i in 1:folds) {
       indX1 = which(indX==i)
       indX2 = which(indX!=i)
       indK1 = which(indKnowns==i)
       indK2 = which(indKnowns!=i)
       
       modelTmp = fun(X = X[indX2,,drop=F], knowns = knowns[indK2,,drop=F], class=class[indK2], k=k, B=B[indK2,,drop=F], P=P[indK2,,drop=F], model.structure=model.structure, ...)
       predTmp  = predict(modelTmp, knowns[indK1,,drop=F])
       if (!is.null(class)) {
           errors[i] = mean(class[indK1,drop=F] != predTmp$class.X)
       } else {
         if (!is.null(B)) {
           errors[i] = mean(abs(B[indK1,,drop=F] - predTmp$tij.X))
         } else {
           errors[i] = mean(abs(P[indK1,,drop=F] - predTmp$tij.X))
         }
       }
   }
   list(errors=errors, indKnowns=indKnowns, indX=indX)
}


