generate.partition <- function(ntot, nfold=5) {
  foldsize <- rep(round(ntot/nfold), nfold-1)
  foldsize <- c(foldsize, ntot-sum(foldsize))
  
  remain <- 1:ntot
  folds <- rep(NA, ntot)
  for (n in 1:(nfold-1)) {
    idxtmp <- sample(remain, size=foldsize[n])
    folds[idxtmp] <- n
    remain <- setdiff(remain, idxtmp)
  }
  folds[remain] <- nfold
  return (as.integer(folds))
}

generate.partitions <- function(npart=1, ntot, nfold=5, ids=1:npart) {
  if (npart>1) ret <- sapply(1:npart, function(x,ntot,nfold) generate.partition(ntot,nfold), ntot, nfold)
  else ret <- as.matrix(generate.partition(ntot,nfold))
  colnames(ret) <- ids[1:npart]
  return (ret)
}

rmse.error <- function(a,b) sqrt(mean((a-b)^2))

# credit (and potential blame!) for this function goes to anonymous stackoverflow user!
load.object <- function(file) {
  env <- new.env()
  load(file, envir=env)
  loadedObjects <- objects(env, all.names=TRUE)
  stopifnot(length(loadedObjects)==1)
  return (env[[loadedObjects]])
}

extract.baselearner.name <- function(config, type="regression") {
  if (type=="regression") return (strsplit(class(config), "[.]Regression.Config")[[1]][1])
  stop("invalid type")
}

regression.extract.response <- function(formula, data) {
  ret <- try(data[,all.vars(formula)[1]], silent=T)
  if (inherits(ret, "try-error")) stop("response extraction failed")
  return (ret)
}

