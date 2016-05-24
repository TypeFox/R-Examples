bsgw.control <- function(scalex=TRUE, iter=1000, burnin=round(iter/2), sd.thresh=1e-4, lambda=0.0, lambdas=lambda, nskip=round(iter/10), alpha.min=0.1, alpha.max=10.0
                         , beta.max=log(20), betas.max=5.0, memlim.gb=8) {
  return (list(scalex=scalex, iter=iter, burnin=burnin, sd.thresh=sd.thresh, lambda=lambda, lambdas=lambdas, nskip=nskip
               , alpha.min=alpha.min, alpha.max=alpha.max, beta.max=beta.max, betas.max=betas.max, memlim.gb=memlim.gb))
}

bsgw.strip.formula <- function(survformula) {
  allvars <- all.vars(survformula)
  return (formula(paste("Surv(", allvars[1], ",", allvars[2], ")~1", sep="")))
}

bsgw.empty.plot <- function(...) {
  plot(0,0,type="l", xlab="", ylab="",...)}

bsgw.calc.pval <- function(x, ref=0.0, na.rm = FALSE) { # add flag for one-sided vs. two-sided
  if (na.rm) x <- x[!is.na(x)]
  bigger <- median(x)>ref
  if (sd(x)<.Machine$double.eps) {
    ret <- NA
  } else {
    ret <- 2*length(which(if (bigger) x<ref else x>ref))/length(x)
  }
  attr(ret, "bigger") <- bigger
  return (ret)
}

bsgw.scale <- function(X, apply.sc, ...) {
  if (missing(apply.sc)) apply.sc <- which(sapply(1:ncol(X), function(n) length(unique(X[,n]))>2))
  ret <- scale(X[,apply.sc], ...)
  X[,apply.sc] <- ret
  attr(X, "centerVec") <- attr(ret, "scaled:center")
  attr(X, "scaleVec") <- attr(ret, "scaled:scale")
  attr(X, "apply.scale") <- apply.sc
  return (X)
}

bsgw.generate.folds <- function(ntot, nfold=5) {
  # determine size of each fold
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
  return (folds)
}

bsgw.generate.folds.eventbalanced <- function(formula, data, nfold=5) {
  statusCol <- all.vars(formula)[2]
  
  index_with_event <- which(data[,statusCol]==1); nwith <- length(index_with_event)
  index_without_event <- which(data[,statusCol]==0); nwithout <- length(index_without_event)
  ret_with_event <- bsgw.generate.folds(nwith, nfold)
  ret_without_event <- bsgw.generate.folds(nwithout, nfold)
  ret_all <- list()
  ret_flat <- rep(NA, nrow(data))
  for (n in 1:nfold) {
    ret_all[[n]] <- c(index_with_event[which(ret_with_event==n)], index_without_event[which(ret_without_event==n)])
    ret_flat[ret_all[[n]]] <- n
  }
  return (ret_flat)
}

