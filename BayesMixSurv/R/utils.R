bayesmixsurv.strip.formula <- function(survformula) {
  allvars <- all.vars(survformula)
  return (formula(paste("Surv(", allvars[1], ",", allvars[2], ")~1", sep="")))
}

bayesmixsurv.control <- function(single=FALSE, alpha2.fixed=NULL, alpha.boundary=1.0, lambda1=1.0, lambda2=lambda1, iter=1000, burnin=round(iter/2)
                                 , sd.thresh=1e-4, scalex=TRUE, nskip=round(iter/10)) {
  return (list(single=single, alpha2.fixed=alpha2.fixed, alpha.boundary=alpha.boundary, lambda1=lambda1, lambda2=lambda2
               , iter=iter, burnin=burnin, sd.thresh=sd.thresh, scalex=scalex, nskip=nskip))
}

bayesmixsurv.scale <- function(X, apply.sc, ...) {
  if (missing(apply.sc)) apply.sc <- which(sapply(1:ncol(X), function(n) length(unique(X[,n]))>2))
  ret <- scale(X[,apply.sc], ...)
  X[,apply.sc] <- ret
  attr(X, "centerVec") <- attr(ret, "scaled:center")
  attr(X, "scaleVec") <- attr(ret, "scaled:scale")
  attr(X, "apply.scale") <- apply.sc
  return (X)
}

bayesmixsurv.calc.pval <- function(x, ref=0.0, na.rm = FALSE) { # add flag for one-sided vs. two-sided
  if (na.rm) x <- x[!is.na(x)]
  bigger <- median(x)>ref
  if (sd(x)<.Machine$double.eps) {
    ret <- NA
  } else {
    ret <- max(2*length(which(if (bigger) x<ref else x>ref))/length(x), 1/length(x)) # max adjustment inspired by MCMCglmm package
  }
  attr(ret, "bigger") <- bigger
  return (ret)
}

bayesmixsurv.empty.plot <- function(...) {
  plot(0,0,type="l", xlab="", ylab="",...)
}

bayesmixsurv.generate.folds <- function(ntot, nfold=5) {
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

bayesmixsurv.generate.folds.eventbalanced <- function(formula, data, nfold=5) {
  statusCol <- all.vars(formula)[2]
  status_levels <- unique(data[,statusCol])
  if (length(status_levels)>2) stop("status field is not binary")
  
  index_with_event <- which(data[,statusCol]==status_levels[1]); nwith <- length(index_with_event)
  index_without_event <- which(data[,statusCol]==status_levels[2]); nwithout <- length(index_without_event)
  ret_with_event <- bayesmixsurv.generate.folds(nwith, nfold)
  ret_without_event <- bayesmixsurv.generate.folds(nwithout, nfold)
  ret_all <- list()
  ret_flat <- rep(NA, nrow(data))
  for (n in 1:nfold) {
    ret_all[[n]] <- c(index_with_event[which(ret_with_event==n)], index_without_event[which(ret_without_event==n)])
    ret_flat[ret_all[[n]]] <- n
  }
  return (ret_flat)
}

