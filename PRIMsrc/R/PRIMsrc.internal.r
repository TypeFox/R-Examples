##########################################################################################################################################
# PRIMsrc
##########################################################################################################################################

##########################################################################################################################################
# SURVIVAL INTERNAL SUBROUTINES
# (never to be called by end-user)
##########################################################################################################################################

##########################################################################################################################################
################
# Usage         :
################
#                   cv.presel(x, times, status,
#                             vs, n, p, K, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.presel <- function(x, times, status,
                      vs, n, p, K, seed) {

  if (vs) {
    digits <- getOption("digits")
    if (is.null(seed)) {
        seed <- floor(runif(n=1, min=0, max=1) * 10^(min(digits,9)))
    } else {
        set.seed(seed)
    }
    nfolds <- max(3,K)
    folds <- cv.folds(y=status, K=nfolds, seed=seed)
    foldid <- as.numeric(folds$which[folds$permkey])
    enalpha <- seq(from=0, to=1, length.out=10)
    lenalpha <- length(enalpha)
    enlambda <- vector(mode="list", length=lenalpha)
    cv.errmu <- vector(mode="list", length=lenalpha)
    cv.errsd <- vector(mode="list", length=lenalpha)
    for (i in 1:lenalpha) {
        cv.fit <- cv.glmnet(x=x, y=Surv(times, status), alpha=enalpha[i], nfolds=nfolds, foldid=foldid, family="cox", maxit=1e5)
        cv.errmu[[i]] <- cv.fit$cvm
        cv.errsd[[i]] <- cv.fit$cvsd
        enlambda[[i]] <- cv.fit$lambda
    }
    cv.errmu <- list2mat(list=cv.errmu, coltrunc="max", fill=NA)
    cv.errsd <- list2mat(list=cv.errsd, coltrunc="max", fill=NA)
    w <- as.numeric(which(x=(cv.errmu == as.numeric(cv.errmu)[which.min(cv.errmu)]), arr.ind=TRUE, useNames=FALSE))
    ww <- as.matrix(which(x=cv.errmu[w[1],w[2]] + cv.errsd[w[1],w[2]] <= cv.errmu - cv.errsd, arr.ind=TRUE, useNames=TRUE))
    wm <- ww - rep.mat(t(w), 2, nrow(ww))
    wm <- w + wm[which.min(apply(abs(wm), 1, sum)),]        # Nearest neighbor to index minimizer
    if (is.empty(wm)) {
        selected <- NULL
        varsign <- NULL
        success <- FALSE
    } else {
        enlambda.min <- enlambda[[wm[1]]][wm[2]]
        enalpha.min <- enalpha[wm[1]]
        fit <- glmnet(x=x, y=Surv(times, status), alpha=enalpha.min, family="cox", maxit=1e5)
        cv.coef <- as.numeric(coef(fit, s=enlambda.min))
        selected <- which(!(is.na(cv.coef)) & (cv.coef != 0))
        if (is.empty(selected)) {
            selected <- NULL
            varsign <- NULL
            success <- FALSE
        } else {
            names(selected) <- colnames(x)[selected]
            cv.coef <- cv.coef[selected]
            varsign <- sign(cv.coef)
            names(varsign) <- colnames(x)[selected]
            success <- TRUE
        }
    }
  } else {
    fit <- glmnet(x=x, y=Surv(times, status), alpha=0, family="cox", maxit=1e5)
    cv.coef <- as.numeric(coef(fit, s=0))
    selected <- which(!(is.na(cv.coef)) & (cv.coef != 0))
    if (is.empty(selected)) {
      selected <- NULL
      varsign <- NULL
      success <- FALSE
    } else {
      names(selected) <- colnames(x)[selected]
      cv.coef <- cv.coef[selected]
      varsign <- sign(cv.coef)
      names(varsign) <- colnames(x)[selected]
      success <- TRUE
    }
  }

  return(list("selected"=selected,
              "varsign"=varsign,
              "success"=success))
}
##########################################################################################################################################





##########################################################################################################################################
################
# Usage         :
################
#                   cv.box.rep(x, times, status,
#                              B, K, arg,
#                              cvtype, decimals,
#                              probval, timeval,
#                              varsign, initcutpts,
#                              parallel, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.box.rep <- function(x, times, status,
                       B, K, arg,
                       cvtype, decimals,
                       probval, timeval,
                       varsign, initcutpts,
                       parallel, seed) {

  CV.maxsteps <- numeric(B)
  CV.boxind <- vector(mode="list", length=B)
  CV.boxcut <- vector(mode="list", length=B)
  CV.support <- vector(mode="list", length=B)
  CV.lhr <- vector(mode="list", length=B)
  CV.lrt <- vector(mode="list", length=B)
  CV.cer <- vector(mode="list", length=B)
  CV.trace <- vector(mode="list", length=B)
  CV.time.bar <- vector(mode="list", length=B)
  CV.prob.bar <- vector(mode="list", length=B)
  CV.max.time.bar <- vector(mode="list", length=B)
  CV.min.prob.bar <- vector(mode="list", length=B)
  b <- 1
  k <- 0
  success <- TRUE
  while (b <= B) {
    cat("replicate : ", b, "\n", sep="")
    if (!parallel) {
      set.seed(seed[b])
      cat("seed : ", seed[b], "\n", sep="")
    }
    if (cvtype == "averaged") {
      CVBOX <- cv.ave.box(x=x, times=times, status=status,
                          K=K, arg=arg, decimals=decimals,
                          probval=probval, timeval=timeval,
                          varsign=varsign, initcutpts=initcutpts,
                          seed=seed[b])
    } else if (cvtype == "combined") {
      CVBOX <- cv.comb.box(x=x, times=times, status=status,
                           K=K, arg=arg, decimals=decimals,
                           probval=probval, timeval=timeval,
                           varsign=varsign, initcutpts=initcutpts,
                           seed=seed[b])
    } else if (cvtype == "none") {
      CVBOX <- cv.comb.box(x=x, times=times, status=status,
                           K=1, arg=arg, decimals=decimals,
                           probval=probval, timeval=timeval,
                           varsign=varsign, initcutpts=initcutpts,
                           seed=seed[b])
    } else {
      stop("Invalid CV type option \n")
    }
    if (!CVBOX$drop) {
      CV.maxsteps[b] <- CVBOX$cvfit$cv.maxsteps
      CV.trace[[b]] <- CVBOX$cvfit$cv.trace
      CV.boxind[[b]] <- CVBOX$cvfit$cv.boxind
      CV.boxcut[[b]] <- CVBOX$cvfit$cv.boxcut
      CV.support[[b]] <- CVBOX$cvfit$cv.stats$cv.support
      CV.lhr[[b]] <- CVBOX$cvfit$cv.stats$cv.lhr
      CV.lrt[[b]] <- CVBOX$cvfit$cv.stats$cv.lrt
      CV.cer[[b]] <- CVBOX$cvfit$cv.stats$cv.cer
      CV.time.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.time.bar
      CV.prob.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.prob.bar
      CV.max.time.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.max.time.bar
      CV.min.prob.bar[[b]] <- CVBOX$cvfit$cv.stats$cv.min.prob.bar
      b <- b + 1
      k <- 0
    } else {
      cat("Could not find one step in at least one of the folds within replicate #", b ,". Retrying replicate with new seed\n", sep="")
      seed[b] <- seed[b] + B
      k <- k + 1
      if (k == B) {
        cat("Could not complete requested replications after ", B ," successive trials. Exiting replications\n", sep="")
        b <- B
        success <- FALSE
      }
    }
  }

  return(list("cv.maxsteps"=CV.maxsteps,
              "cv.trace"=CV.trace,
              "cv.boxind"=CV.boxind,
              "cv.boxcut"=CV.boxcut,
              "cv.support"=CV.support,
              "cv.lhr"=CV.lhr,
              "cv.lrt"=CV.lrt,
              "cv.cer"=CV.cer,
              "cv.time.bar"=CV.time.bar,
              "cv.prob.bar"=CV.prob.bar,
              "cv.max.time.bar"=CV.max.time.bar,
              "cv.min.prob.bar"=CV.min.prob.bar,
              "success"=success,
              "seed"=seed))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    cv.pval (x, times, status,
#                             cvtype, decimals,
#                             varsign, initcutpts,
#                             A, K, arg, obs.chisq,
#                             parallel, conf)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.pval <- function(x, times, status,
                    cvtype, decimals,
                    varsign, initcutpts,
                    A, K, arg, obs.chisq,
                    parallel, conf) {

  if (!parallel) {
    null.chisq <- cv.null(x=x, times=times, status=status,
                          cvtype=cvtype, decimals=decimals,
                          varsign=varsign, initcutpts=initcutpts,
                          A=A, K=K, arg=arg)
  } else {
    if (conf$type == "SOCK") {
      cl <- makeCluster(spec=conf$names,
                        type=conf$type,
                        homogeneous=conf$homo,
                        outfile=conf$outfile,
                        verbose=conf$verbose)
    } else {
      cl <- makeCluster(spec=conf$cpus,
                        type=conf$type,
                        homogeneous=conf$homo,
                        outfile=conf$outfile,
                        verbose=conf$verbose)
    }
    clusterSetRNGStream(cl=cl, iseed=NULL)
    null.cl <- clusterCall(cl=cl, fun=cv.null,
                           x=x, times=times, status=status,
                           cvtype=cvtype, decimals=decimals,
                           varsign=varsign, initcutpts=initcutpts,
                           A=ceiling(A/conf$cpus), K=K, arg=arg)
    stopCluster(cl)
    null.chisq <- cbindlist(null.cl)
  }
  cvl <- nrow(null.chisq)
  pval <- numeric(cvl)
  for (l in 1:cvl) {
    pval[l] <- mean((null.chisq[l,] >= obs.chisq[l]), na.rm=TRUE)
  }
  pval <- round(pval, digits=floor(log(base=10, A)))
  names(pval) <- paste("step", 0:(cvl-1), sep="")

  return(pval)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                   cv.null (x, times, status,
#                            cvtype, decimals,
#                            varsign, initcutpts,
#                            A, K, arg)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.null <- function(x, times, status,
                    cvtype, decimals,
                    varsign, initcutpts,
                    A, K, arg) {

  n <- nrow(x)
  null.chisq <- vector(mode="list", length=A)
  a <- 1
  while (a <= A) {
    cat("Permutation sample: ", a, "\n")
    perm.ind <- sample(x = 1:n, size = n, replace = FALSE, prob = NULL)
    perm.times <- times[perm.ind]
    perm.status <- status[perm.ind]
    if (cvtype == "averaged") {
      obj <- tryCatch({cv.ave.box(x=x, times=perm.times, status=perm.status,
                                  varsign=varsign, initcutpts=initcutpts,
                                  K=K, arg=arg, decimals=decimals,
                                  probval=NULL, timeval=NULL, seed=NULL)}, error=function(w){NULL})
      if (is.list(obj)) {
        null.chisq[[a]] <- obj$cvfit$cv.stats$cv.lrt
        a <- a + 1
      } else {
        cat("Permutation sample dropped... \n")
      }
    } else if (cvtype == "combined") {
      obj <- tryCatch({cv.comb.box(x=x, times=perm.times, status=perm.status,
                                   varsign=varsign, initcutpts=initcutpts,
                                   K=K, arg=arg, decimals=decimals,
                                   probval=NULL, timeval=NULL, seed=NULL)}, error=function(w){NULL})
      if (is.list(obj)) {
        null.chisq[[a]] <- obj$cvfit$cv.stats$cv.lrt
        a <- a + 1
      } else {
        cat("Permutation sample dropped... \n")
      }
    } else if (cvtype == "none") {
      obj <- tryCatch({cv.comb.box(x=x, times=perm.times, status=perm.status,
                                   varsign=varsign, initcutpts=initcutpts,
                                   K=1, arg=arg, decimals=decimals,
                                   probval=NULL, timeval=NULL, seed=NULL)}, error=function(w){NULL})
      if (is.list(obj)) {
        null.chisq[[a]] <- obj$cvfit$cv.stats$cv.lrt
        a <- a + 1
      } else {
        cat("Permutation sample dropped... \n")
      }
    } else {
      stop("Invalid CV type option \n")
    }
  }

  return(t(list2mat(null.chisq)))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                   cv.ave.box (x, times, status,
#                               probval, timeval, decimals,
#                               varsign, initcutpts,
#                               K, arg, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.ave.box <- function(x, times, status,
                       probval, timeval, decimals,
                       varsign, initcutpts,
                       K, arg, seed) {

  n <- nrow(x)
  p <- ncol(x)

  alpha <- NULL
  beta <- NULL
  minn <- NULL
  L <- NULL
  peelcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=arg, split=",")) ))

  fold.obj <- cv.ave.fold(x=x, times=times, status=status,
                          probval=probval, timeval=timeval,
                          varsign=varsign, initcutpts=initcutpts,
                          K=K, arg=arg, seed=seed)
  trace.list <- fold.obj$trace.list
  boxstat.list <- fold.obj$boxstat.list
  boxcut.list <- fold.obj$boxcut.list
  drop <- fold.obj$drop

  # Cross-validated maximum peeling length from all folds
  CV.Lm <- min(fold.obj$nsteps)

  # Truncate the cross-validated quantities from all folds to the same cross-validated length
  for (k in 1:K) {
    boxcut.list[[k]] <- boxcut.list[[k]][1:CV.Lm,,drop=FALSE]
  }

  # Compute the averaged box statistics for each step from all the folds
  # Each entry or row signifies a step
  CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(x)))
  for (l in 1:CV.Lm) {
    summincut <- matrix(NA, K, p)
    for (k in 1:K) {
        summincut[k,] <- boxcut.list[[k]][l,]
    }
    CV.boxcut[l, ] <- colMeans(summincut, na.rm=TRUE)
  }
  rownames(CV.boxcut) <- paste("step", 0:(CV.Lm-1), sep="")
  colnames(CV.boxcut) <- colnames(x)
  
  # Get the box membership indicator vector of all observations for each step from all the folds
  # Based on the corresponding averaged box over the folds
  CV.boxind <- matrix(NA, nrow=CV.Lm, ncol=n)
  for (l in 1:CV.Lm) {
    boxcut <- CV.boxcut[l, ] * varsign
    x.cut <- t(t(x) * varsign)
    x.ind <- t(t(x.cut) >= boxcut)
    CV.boxind[l,] <- (rowMeans(x.ind) == 1)  # Set as TRUE which observations are inside the box boudaries for all axes directions
  }
  rownames(CV.boxind) <- paste("step", 0:(CV.Lm-1), sep="")
  colnames(CV.boxind) <- rownames(x)

  # Get the adjusted cross-validated maximum peeling length, thresholded by minimal box support
  CV.Lm <- max(which(apply(CV.boxind, 1, function(x) {length(which(x))/n >= max(minn/n, beta)})))
        
  # Get the adjusted test box defnition and membership indicator vector of all observations for each step from all the folds
  CV.boxind <- CV.boxind[1:CV.Lm,,drop=FALSE]
  CV.boxcut <- CV.boxcut[1:CV.Lm,,drop=FALSE]

  # Get the covariates traces from all folds
  # Variable traces are first stacked and truncated in a matrix where folds are by rows and steps by columns
  CV.trace <- list2mat(list=trace.list, coltrunc=CV.Lm)
  CV.trace <- t(CV.trace)
  dimnames(CV.trace) <- list(paste("step", 0:(CV.Lm-1), sep=""), 1:K)

  # Truncate the cross-validated quantities from all folds to the same cross-validated length
  for (k in 1:K) {
    boxstat.list[[k]] <- boxstat.list[[k]][1:CV.Lm]
  }
  
  # Compute the averaged box statistics for each step from all the folds
  # Each entry or row signifies a step
  CV.support <- rep(NA, CV.Lm)
  names(CV.support) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.lhr <- rep(NA, CV.Lm)
  names(CV.lhr) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.lrt <- rep(NA, CV.Lm)
  names(CV.lrt) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.cer <- rep(NA, CV.Lm)
  names(CV.cer) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.time.bar <- rep(NA, CV.Lm)
  names(CV.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.prob.bar <- rep(NA, CV.Lm)
  names(CV.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.max.time.bar <- rep(NA, CV.Lm)
  names(CV.max.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.min.prob.bar <- rep(NA, CV.Lm)
  names(CV.min.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  for (l in 1:CV.Lm) {
    sumtime <- rep(NA, K)
    sumprob <- rep(NA, K)
    summaxtime <- rep(NA, K)
    summinprob <- rep(NA, K)
    sumlhr <- rep(NA, K)
    sumlrt <- rep(NA, K)
    sumcer <- rep(NA, K)
    sumsupport <- rep(NA, K)
    for (k in 1:K) {
      outbounds <- boxstat.list[[k]][[l]]
      if (!is.null(outbounds)) {
        sumlhr[k] <- boxstat.list[[k]][[l]][[1]]
        sumlrt[k] <- boxstat.list[[k]][[l]][[2]]
        sumcer[k] <- boxstat.list[[k]][[l]][[3]]
        sumsupport[k] <- boxstat.list[[k]][[l]][[4]]
        sumtime[k] <- boxstat.list[[k]][[l]][[5]]
        sumprob[k] <- boxstat.list[[k]][[l]][[6]]
        summaxtime[k] <- boxstat.list[[k]][[l]][[7]]
        summinprob[k] <- boxstat.list[[k]][[l]][[8]]
      }
    }
    CV.lhr[l] <- mean(sumlhr, na.rm=TRUE)
    CV.lrt[l] <- mean(sumlrt, na.rm=TRUE)
    CV.cer[l] <- mean(sumcer, na.rm=TRUE)
    CV.support[l] <- mean(sumsupport, na.rm=TRUE)
    CV.time.bar[l] <- mean(sumtime, na.rm=TRUE)
    CV.prob.bar[l] <- mean(sumprob, na.rm=TRUE)
    CV.max.time.bar[l] <- mean(summaxtime, na.rm=TRUE)
    CV.min.prob.bar[l] <- mean(summinprob, na.rm=TRUE)
  }

  # Box peeling rules for each step
  CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(x))))
  for (j in 1:p) {
    if (varsign[j] > 0) {
      ss <- ">="
    } else {
      ss <- "<="
    }
    CV.rules[, j] <- paste(colnames(x)[j], ss, format(x=CV.boxcut[, j], digits=decimals, nsmall=decimals), sep="")
  }

  # Applying the cross-validation criterion to the profiles
  # Cross-validated optimal length from all folds
  # By maximization of the LHR (between in and out box test samples)
  if (all(is.na(CV.lhr))) {
    CV.L.lhr <- NA
  } else {
    CV.L.lhr <- which.max(CV.lhr)
  }
  # By maximization of the LRT (between in and out box test samples)
  if (all(is.na(CV.lrt))) {
    CV.L.lrt <- NA
  } else {
    CV.L.lrt <- which.max(CV.lrt)
  }
  # By minimization of the CER (between predicted and observed inbox test samples survival times)
  if (all(is.na(CV.cer))) {
    CV.L.cer <- NA
  } else {
    CV.L.cer <- which.min(CV.cer)
  }

  # Box statistics for each step
  CV.stats <-  data.frame("cv.support"=CV.support,
                          "cv.lhr"=CV.lhr,
                          "cv.lrt"=CV.lrt,
                          "cv.cer"=CV.cer,
                          "cv.time.bar"=CV.time.bar,
                          "cv.prob.bar"=CV.prob.bar,
                          "cv.max.time.bar"=CV.max.time.bar,
                          "cv.min.prob.bar"=CV.min.prob.bar)
  rownames(CV.stats) <- paste("step", 0:(CV.Lm-1), sep="")

  # Create the return object 'CV.fit'
  CV.fit <- list("cv.nsteps.lhr"=CV.L.lhr,
                 "cv.nsteps.lrt"=CV.L.lrt,
                 "cv.nsteps.cer"=CV.L.cer,
                 "cv.maxsteps"=CV.Lm,
                 "cv.boxcut"=CV.boxcut,
                 "cv.rules"=CV.rules,
                 "cv.stats"=CV.stats,
                 "cv.trace"=CV.trace,
                 "cv.boxind"=CV.boxind)

  return(list("x"=x, "times"=times, "status"=status,
              "cvfit"=CV.fit, "drop"=drop, "seed"=seed))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                   cv.comb.box (x, times, status,
#                                probval, timeval, decimals,
#                                varsign, initcutpts,
#                                K, arg, seed)
#
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.comb.box <- function(x, times, status,
                        probval, timeval, decimals,
                        varsign, initcutpts,
                        K, arg, seed) {
  n <- nrow(x)
  p <- ncol(x)

  alpha <- NULL
  beta <- NULL
  minn <- NULL
  L <- NULL
  peelcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=arg, split=",")) ))

  fold.obj <- cv.comb.fold(x=x, times=times, status=status,
                           varsign=varsign, initcutpts=initcutpts,
                           K=K, arg=arg, seed=seed)
  
  ord <- fold.obj$key
  times.list <- fold.obj$cvtimes
  status.list <- fold.obj$cvstatus
  trace.list <- fold.obj$trace
  boxind.list <- fold.obj$boxind
  boxcut.list <- fold.obj$boxcut

  # Cross-validated maximum peeling length from all folds
  CV.Lm <- min(fold.obj$nsteps)

  # Get the test box membership indicator vector of all observations for each step from all the folds
  # Based on the combined membership indicator vectors over the folds
  # Re-ordered by initial order of observations
  CV.boxind <- cbindlist(boxind.list, trunc=CV.Lm)[,ord,drop=FALSE]
  rownames(CV.boxind) <- paste("step", 0:(CV.Lm-1), sep="")
  colnames(CV.boxind) <- rownames(x)

  # Get the adjusted cross-validated maximum peeling length, thresholded by minimal box support
  CV.Lm <- max(which(apply(CV.boxind, 1, function(x) {length(which(x))/n >= max(minn/n, beta)})))
        
  # Get the adjusted test box membership indicator vector of all observations for each step from all the folds
  CV.boxind <- CV.boxind[1:CV.Lm,,drop=FALSE]

  # Concatenates the observations of test times and status from all folds
  # Re-ordered by initial order of observations
  CV.times <- unlist(times.list)[ord]
  CV.status <- unlist(status.list)[ord]

  # Get the variable traces
  # Variable traces are first stacked and truncated in a matrix where folds are by rows and steps by columns
  CV.trace <- list2mat(list=trace.list, coltrunc=CV.Lm)
  CV.trace <- t(CV.trace)
  dimnames(CV.trace) <- list(paste("step", 0:(CV.Lm-1), sep=""), 1:K)

  # Get the combined boxcut (truncated to the same cross-validated length) for each step from all the folds
  # using the circumscribing box to the conmbined test set in-box samples over all the folds
  CV.boxcut <- matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(x)))
  tmparray <- list2array(list=boxcut.list, rowtrunc=CV.Lm)
  for (l in 1:CV.Lm) {
    for (j in 1:p) {
        if (varsign[j] > 0) {
          CV.boxcut[l,j] <- min(x[CV.boxind[l,],j])
        } else {
          CV.boxcut[l,j] <- max(x[CV.boxind[l,],j])
        }
    }
  }

  # Box peeling rules for each step
  CV.rules<- as.data.frame(matrix(data=NA, nrow=CV.Lm, ncol=p, dimnames=list(paste("step", 0:(CV.Lm-1), sep=""), colnames(x))))
  for (j in 1:p) {
    if (varsign[j] > 0) {
      ss <- ">="
    } else {
      ss <- "<="
    }
    CV.rules[, j] <- paste(colnames(x)[j], ss, format(x=CV.boxcut[, j], digits=decimals, nsmall=decimals), sep="")
  }

  # Compute the combined test box statistics from all folds for all steps, each entry or row signifies a step
  CV.support <- rep(NA, CV.Lm)
  names(CV.support) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.lhr <- rep(NA, CV.Lm)
  names(CV.lhr) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.lrt <- rep(NA, CV.Lm)
  names(CV.lrt) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.cer <- rep(NA, CV.Lm)
  names(CV.cer) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.time.bar <- rep(NA, CV.Lm)
  names(CV.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.prob.bar <- rep(NA, CV.Lm)
  names(CV.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.max.time.bar <- rep(NA, CV.Lm)
  names(CV.max.time.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  CV.min.prob.bar <- rep(NA, CV.Lm)
  names(CV.min.prob.bar) <- paste("step", 0:(CV.Lm-1), sep="")
  timemat <- matrix(NA, nrow=CV.Lm, ncol=n)
  probmat <- matrix(NA, nrow=CV.Lm, ncol=n)
  ind.rem <- numeric(0)
  for (l in 1:CV.Lm) {
    boxind <- CV.boxind[l,]
    boxind1 <- 1*boxind
    if ((l == 1) && (sum(boxind, na.rm=TRUE) != 0)) {
      surv.fit <- survfit(Surv(CV.times[boxind], CV.status[boxind]) ~ 1)
      timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
      probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
      CV.lhr[l] <- 0
      CV.lrt[l] <- 0
      CV.cer[l] <- 1
      CV.support[l] <- 1
    } else if ((sum(boxind, na.rm=TRUE) != length(boxind[!is.na(boxind)])) && (sum(boxind, na.rm=TRUE) != 0)) {
      surv.fit <- survfit(Surv(CV.times[boxind], CV.status[boxind]) ~ 1)
      timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
      probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
      surv.formula <- (Surv(CV.times, CV.status) ~ 1 + boxind1)
      coxobj <- coxph(surv.formula, singular.ok=TRUE, iter.max=1)
      CV.lhr[l] <- coxobj$coef
      CV.lrt[l] <- survdiff(surv.formula, rho=0)$chisq
      predobj <- predict(object=coxobj, type="lp", reference="sample")
      CV.cer[l] <- rcorr.cens(x=predobj, S=Surv(CV.times, CV.status))['C Index']
      CV.support[l] <- mean(boxind, na.rm=TRUE)
    } else {
      timemat[l, ] <- NA
      probmat[l, ] <- NA
      CV.lhr[l] <- 0
      CV.lrt[l] <- 0
      CV.cer[l] <- 1
      CV.support[l] <- NA
      ind.rem <- c(ind.rem, l)
    }
  }
  if (length(ind.rem) != CV.Lm) {
    drop <- FALSE
    endobj <- endpoints (ind=ind.rem, timemat=timemat, probmat=probmat, timeval=timeval, probval=probval)
    time.bar <- endobj$time.bar
    prob.bar <- endobj$prob.bar
    max.time.bar <- endobj$max.time.bar
    min.prob.bar <- endobj$min.prob.bar
    for (l in 1:CV.Lm) {
      if (!(l %in% ind.rem)) {
        CV.time.bar[l] <- time.bar[l]
        CV.prob.bar[l] <- prob.bar[l]
        CV.max.time.bar[l] <-  max.time.bar[l]
        CV.min.prob.bar[l] <- min.prob.bar[l]
      } else {
        CV.time.bar[l] <- NA
        CV.prob.bar[l] <- NA
        CV.max.time.bar[l] <- NA
        CV.min.prob.bar[l] <- NA
      }
    }
  } else {
    cat("Dropped !\n", sep="")
    drop <- TRUE
    CV.time.bar <- rep(NA, CV.Lm)
    CV.prob.bar <- rep(NA, CV.Lm)
    CV.max.time.bar <- rep(NA, CV.Lm)
    CV.min.prob.bar <- rep(NA, CV.Lm)
  }

  # Applying the cross-validation criterion to the profiles
  # Cross-validated optimal length from all folds
  # By maximization of the LHR (between in and out box test samples)
  if (all(is.na(CV.lhr))) {
    CV.L.lhr <- NA
  } else {
    CV.L.lhr <- which.max(CV.lhr)
  }
  # By maximization of the LRT (between in and out box test samples)
  if (all(is.na(CV.lrt))) {
    CV.L.lrt <- NA
  } else {
    CV.L.lrt <- which.max(CV.lrt)
  }
  # By minimization of the CER (between predicted and observed inbox test samples survival times)
  if (all(is.na(CV.cer))) {
    CV.L.cer <- NA
  } else {
    CV.L.cer <- which.min(CV.cer)
  }

  # Box statistics for each step
  CV.stats <-  data.frame("cv.support"=CV.support,
                          "cv.lhr"=CV.lhr,
                          "cv.lrt"=CV.lrt,
                          "cv.cer"=CV.cer,
                          "cv.time.bar"=CV.time.bar,
                          "cv.prob.bar"=CV.prob.bar,
                          "cv.max.time.bar"=CV.max.time.bar,
                          "cv.min.prob.bar"=CV.min.prob.bar)
  rownames(CV.stats) <- paste("step", 0:(CV.Lm-1), sep="")

  # Create the return object 'CV.fit'
  CV.fit <- list("cv.nsteps.lhr"=CV.L.lhr,
                 "cv.nsteps.lrt"=CV.L.lrt,
                 "cv.nsteps.cer"=CV.L.cer,
                 "cv.maxsteps"=CV.Lm,
                 "cv.boxcut"=CV.boxcut,
                 "cv.rules"=CV.rules,
                 "cv.stats"=CV.stats,
                 "cv.trace"=CV.trace,
                 "cv.boxind"=CV.boxind)

  return(list("x"=x, "times"=times, "status"=status,
              "cvfit"=CV.fit, "drop"=drop, "seed"=seed))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    cv.ave.fold (x, times, status,
#                                 probval, timeval,
#                                 varsign, initcutpts,
#                                 K, arg, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.ave.fold <- function(x, times, status,
                        probval, timeval,
                        varsign, initcutpts,
                        K, arg, seed) {

  drop <- FALSE
  folds <- cv.folds(y=status, K=K, seed=seed)

  boxstat.list <- vector(mode="list", length=K)
  boxcut.list <- vector(mode="list", length=K)
  trace.list <- vector(mode="list", length=K)
  nsteps <- numeric(K)

  for (k in 1:K) {
    cat("Fold : ", k, "\n", sep="")
    # Initialize training and test data
    if (K == 1) {
      traindata <- testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- times[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- status[folds$perm[(folds$which == k)]]
    } else {
      traindata <- x[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- times[folds$perm[(folds$which != k)]]
      trainstatus <- status[folds$perm[(folds$which != k)]]
      testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- times[folds$perm[(folds$which == k)]]
      teststatus <- status[folds$perm[(folds$which == k)]]
    }
    peelobj <- cv.ave.peel(traindata=traindata, trainstatus=trainstatus, traintime=traintime,
                           testdata=testdata, teststatus=teststatus, testtime=testtime,
                           probval=probval, timeval=timeval,
                           varsign=varsign, initcutpts=initcutpts, arg=arg)
    # Store the test set results from each fold
    nsteps[k] <- peelobj$nsteps
    boxstat.list[[k]] <- peelobj$boxstat
    boxcut.list[[k]] <- peelobj$boxcut
    trace.list[[k]] <- peelobj$trace
    drop <- (drop || peelobj$drop)
  }

  return(list("nsteps"=nsteps,
              "boxstat.list"=boxstat.list,
              "boxcut.list"=boxcut.list,
              "trace.list"=trace.list,
              "drop"=drop))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    cv.comb.fold (x, times, status,
#                                  varsign, initcutpts,
#                                  K, arg, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.comb.fold <- function(x, times, status,
                         varsign, initcutpts,
                         K, arg, seed) {

  folds <- cv.folds(y=status, K=K, seed=seed)
  cvtimes <- vector(mode="list", length=K)
  cvstatus <- vector(mode="list", length=K)
  boxind <- vector(mode="list", length=K)
  boxcut <- vector(mode="list", length=K)
  trace <- vector(mode="list", length=K)
  nsteps <- numeric(K)

  for (k in 1:K) {
    cat("Fold : ", k, "\n", sep="")
    if (K == 1) {
      traindata <- testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      traintime <- testtime <- times[folds$perm[(folds$which == k)]]
      trainstatus <- teststatus <- status[folds$perm[(folds$which == k)]]
    } else {
      traindata <- x[folds$perm[(folds$which != k)], , drop=FALSE]
      traintime <- times[folds$perm[(folds$which != k)]]
      trainstatus <- status[folds$perm[(folds$which != k)]]
      testdata <- x[folds$perm[(folds$which == k)], , drop=FALSE]
      testtime <- times[folds$perm[(folds$which == k)]]
      teststatus <- status[folds$perm[(folds$which == k)]]
    }
    # Store the test set data from each fold (Note: the order of observations of times and status from each fold is kept in the list)
    cvtimes[[k]] <- testtime
    cvstatus[[k]] <- teststatus
    peelobj <- cv.comb.peel(traindata=traindata, trainstatus=trainstatus, traintime=traintime,
                            testdata=testdata, teststatus=teststatus, testtime=testtime,
                            varsign=varsign, initcutpts=initcutpts, arg=arg)
    # Store the test set results from each fold
    nsteps[k] <- peelobj$nsteps
    boxind[[k]] <- peelobj$boxind
    boxcut[[k]] <- peelobj$boxcut
    trace[[k]] <- peelobj$trace
  }

  return(list("nsteps"=nsteps, "key"=folds$foldkey,
              "cvtimes"=cvtimes, "cvstatus"=cvstatus,
              "boxind"=boxind, "boxcut"=boxcut, "trace"=trace))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    cv.ave.peel (traindata, trainstatus, traintime,
#                                 testdata, teststatus, testtime,
#                                 probval, timeval,
#                                 varsign, initcutpts, arg)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.ave.peel <- function(traindata, trainstatus, traintime,
                        testdata, teststatus, testtime,
                        probval, timeval,
                        varsign, initcutpts, arg) {

  # Training the model
  peelobj <- peel.box(traindata=traindata, traintime=traintime, trainstatus=trainstatus,
                      varsign=varsign, initcutpts=initcutpts, arg=arg)
  nsteps <- peelobj$nsteps

  # Compute the box statistics for all steps, each entry or row signifies a step
  boxstat <- vector(mode="list", length=nsteps)
  timemat <- matrix(NA, nrow=nsteps, ncol=nrow(testdata))
  probmat <- matrix(NA, nrow=nsteps, ncol=nrow(testdata))
  ind.rem <- numeric(0)
  for (l in 1:nsteps) {
    # Extract the rule and sign as one vector
    boxcut <- peelobj$boxcut[l, ] * varsign
    test.cut <- t(t(testdata) * varsign)
    test.ind <- t(t(test.cut) >= boxcut)
    # Set as TRUE which observations are TRUE for all covariates
    test.ind <- (rowMeans(test.ind) == 1)
    test.ind1 <- 1*test.ind
    if ((l == 1) && (sum(test.ind, na.rm=TRUE) != 0)) {
      lhr <- 0
      lrt <- 0
      cer <- 1
      support <- 1
      boxstat[[l]] <- c(lhr, lrt, cer, support)
      names(boxstat[[l]]) <- NULL
      surv.fit <- survfit(Surv(testtime[test.ind], teststatus[test.ind]) ~ 1)
      timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
      probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
    } else if ((sum(test.ind, na.rm=TRUE) != length(test.ind[!is.na(test.ind)])) && (sum(test.ind, na.rm=TRUE) != 0)) {
      surv.formula <- (Surv(testtime, teststatus) ~ 1 + test.ind1)
      coxobj <- coxph(surv.formula, singular.ok=TRUE, iter.max=1)
      lhr <- coxobj$coef
      lrt <- survdiff(surv.formula, rho=0)$chisq
      predobj <- predict(object=coxobj, type="lp", reference="sample")
      cer <- rcorr.cens(x=predobj, S=Surv(testtime, teststatus))['C Index']
      support <- mean(test.ind, na.rm=TRUE)
      boxstat[[l]] <- c(lhr, lrt, cer, support)
      names(boxstat[[l]]) <- NULL
      surv.fit <- survfit(Surv(testtime[test.ind], teststatus[test.ind]) ~ 1)
      timemat[l, (1:length(surv.fit$time))] <- surv.fit$time
      probmat[l, (1:length(surv.fit$surv))] <- surv.fit$surv
    } else {
      lhr <- 0
      lrt <- 0
      cer <- 1
      support <- NA
      boxstat[[l]] <- c(lhr, lrt, cer, support)
      names(boxstat[[l]]) <- NULL
      timemat[l, ] <- NA
      probmat[l, ] <- NA
      ind.rem <- c(ind.rem, l)
    }
  }
  if (length(ind.rem) != nsteps) {
    drop <- FALSE
    endobj <- endpoints (ind=ind.rem, timemat=timemat, probmat=probmat, timeval=timeval, probval=probval)
    time.bar <- endobj$time.bar
    prob.bar <- endobj$prob.bar
    max.time.bar <- endobj$max.time.bar
    min.prob.bar <- endobj$min.prob.bar
    for (l in 1:nsteps) {
      if (!(l %in% ind.rem)) {
        boxstat[[l]] <- c(boxstat[[l]], time.bar[l], prob.bar[l], max.time.bar[l], min.prob.bar[l])
      } else {
        boxstat[[l]] <- c(boxstat[[l]], NA, NA, NA, NA)
      }
    }
  } else {
    cat("Dropped !\n", sep="")
    drop <- TRUE
    for (l in 1:nsteps) {
      boxstat[[l]] <- rep(NA, 8)
    }
  }

  return(list("nsteps"=nsteps, "boxstat"=boxstat, "boxcut"=peelobj$boxcut, "trace"=peelobj$trace, "drop"=drop))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    cv.comb.peel (traindata, trainstatus, traintime,
#                                  testdata, teststatus, testtime,
#                                  varsign, initcutpts, arg)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.comb.peel <- function(traindata, trainstatus, traintime,
                         testdata, teststatus, testtime,
                         varsign, initcutpts, arg) {
  # Training the model
  peelobj <- peel.box(traindata=traindata, traintime=traintime, trainstatus=trainstatus,
                      varsign=varsign, initcutpts=initcutpts, arg=arg)
  nsteps <- peelobj$nsteps

  # Create the indicator matrix of the test data that is within the box for each step
  boxind <- matrix(NA, nrow=nsteps, ncol=nrow(testdata))
  for (l in 1:nsteps) {
    # Extract the rule and sign as one vector
    boxcut <- peelobj$boxcut[l, ] * varsign
    test.cut <- t(t(testdata) * varsign)
    test.ind <- t(t(test.cut) >= boxcut)
    # Set as TRUE which observations are TRUE for all covariates
    boxind[l, ] <- (rowMeans(test.ind) == 1)
  }

  return(list("boxind"=boxind, "nsteps"=nsteps, "boxcut"=peelobj$boxcut, "trace"=peelobj$trace))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    peel.box (traindata, traintime, trainstatus,
#                              varsign, initcutpts, arg)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

peel.box <- function(traindata, traintime, trainstatus,
                     varsign, initcutpts, arg) {

  alpha <- NULL
  beta <- NULL
  minn <- NULL
  L <- NULL
  peelcriterion <- NULL
  eval(parse( text=unlist(strsplit(x=arg, split=",")) ))

  digits <- getOption("digits")

  # Ensures that the training data is a numeric matrix
  traindata <- as.matrix(traindata)
  mode(traindata) <- "numeric"

  # Constants
  n <- nrow(traindata)                                   # Number of samples
  p <- ncol(traindata)                                   # Number of initially pre-selected covariates
  beta <- max(minn/n, beta)                              # Threshold of minimal box support by `minn` points
  ncut <- ceiling(log(1/n) / log(1 - (1/n)))             # Maximal possible number of peeling steps

  # Initializations of variable trace and box boundaries
  trace <- numeric(ncut)
  boxcut <- matrix(data=NA, nrow=ncut, ncol=p)

  # Initializations
  boxcutpts <- initcutpts
  boxmass <- 1
  boxes <- matrix(data=FALSE, nrow=n, ncol=p)            # Initial logical matrix of box membership indicator by dimension
  sel <- rep(TRUE, n)                                    # Initial logical vector of in-box samples
  xsel <- traindata                                      # Initial selection of samples from training data
  varpeel <- (apply(traindata, 2, "var") > 10^(-digits)) # Initial selection of covariates for peeling
  continue <- TRUE
  if (!(is.null(L))) {
    switch <- 1
  } else {
    L <- 1
    switch <- 0
  }
  l <- 0
  lhrlj <- matrix(NA, ncut, p)
  lrtlj <- matrix(NA, ncut, p)
  chslj <- matrix(NA, ncut, p)

  while ((boxmass >= beta) & (l*switch < L) & (continue)) {
    l <- l + 1
    xsign <- t(t(xsel) * varsign)

    # Potential cutpts by dimension
    cutpts.sign <- updatecut(x=xsign, fract=alpha)
    cutpts <- cutpts.sign * varsign

    # Update box membership indicator by dimension
    boxes <- as.matrix(t((t(traindata) * varsign) >= as.vector(cutpts.sign)) & sel)

    vmd <- rep(NA, p)
    for (j in 1:p) {
      boxes1j <- 1 * boxes[,j]
      if ((sum(boxes1j) != length(boxes1j)) && (sum(boxes1j) != 0)) {
        # Rate of increase of LHR (between in and out box)
        if (peelcriterion == "hr") {
          lhrlj[l,j] <- coxph(formula=Surv(traintime, trainstatus) ~ 1 + boxes1j, singular.ok=TRUE, iter.max=1)$coef
          if (l == 1) {
            vmd[j] <- (lhrlj[l,j] - 0) / (1 - mean(boxes1j))
          } else {
            vmd[j] <- (lhrlj[l,j] - lhrlj[l-1,j]) / (boxmass - mean(boxes1j))
          }
        # Rate of increase of LRT (between in and out box)
        } else if (peelcriterion == "lr") {
          lrtlj[l,j] <- survdiff(formula=Surv(traintime, trainstatus) ~ 1 + boxes1j, rho=0)$chisq
          if (l == 1) {
            vmd[j] <- (lrtlj[l,j] - 0) / (1 - mean(boxes1j))
          } else {
            vmd[j] <- (lrtlj[l,j] - lrtlj[l-1,j]) / (boxmass - mean(boxes1j))
          }
        # Rate of increase of CHS (between in and out box)
        } else if (peelcriterion == "ch") {
          fit <- survfit(formula=Surv(traintime, trainstatus) ~ 1, subset=(boxes1j == 1))
          chslj[l,j] <- sum(cumsum(fit$n.event/fit$n.risk))
          if (l == 1) {
            vmd[j] <- (chslj[l,j] - 0) / (1 - mean(boxes1j))
          } else {
            vmd[j] <- (chslj[l,j] - chslj[l-1,j]) / (boxmass - mean(boxes1j))
          }
        } else {
          stop("Invalid peeling criterion \n")
        }
      } else {
        varpeel[j] <- FALSE
      }
    }

    # If the previous attempted peeling succeeded
    if (sum(varpeel) > 0 && (!is.empty(vmd[(!is.nan(vmd)) & (!is.infinite(vmd)) & (!is.na(vmd))]))) {
      # Maximizing the rate of increase of peeling criterion.
      # Only one variable (the first one in rank) is selected in case of ties
      varj <- which(vmd == max(vmd[(!is.nan(vmd)) & (!is.infinite(vmd))], na.rm=TRUE))[1]
      # Updating
      sel <- boxes[, varj, drop=TRUE]
      boxmass <- mean(1 * sel)
      xsel <- traindata[sel, ,drop=FALSE]
      varpeel <- (apply(xsel, 2, "var") > 10^(-digits))
      boxcutpts[varj] <- cutpts[varj]
      # Saving trained box quantities of interest for the current peeling step
      boxcut[l, ] <- boxcutpts
      trace[l] <- varj
    # Else exit the loop and decrement the peeling step number since the last attempted failed in that case
    } else {
      continue <- FALSE
      l <- l - 1
    }
  }

  if (l == 0) {
    # Taking the first step box covering all the data
    boxcut <- rbind(initcutpts)
    trace <- 0
  } else if (l >= 1) {
    # Prepending the first step box covering all the data
    boxcut <- rbind(initcutpts, boxcut[1:l, , drop=FALSE])
    trace <- c(0, trace[1:l])
  }

  rownames(boxcut) <- paste("step", 0:l, sep="")
  colnames(boxcut) <- colnames(traindata)
  names(trace) <- paste("step", 0:l, sep="")

  # Returning the final results, considering the starting point as step #0
  return(list("nsteps"=l+1,
              "boxcut"=boxcut,
              "trace"=trace))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    cv.folds (y, K, seed)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cv.folds <- function (y, K, seed) {

    if (!is.null(seed))
        set.seed(seed)

    n <- length(y)
    K <- round(rep(K, length.out = 1))
    if (!isTRUE((K >= 1) && K <= n))
        stop(paste("`K` is outside the allowable range {1,...,", n, "} \n", sep=""))

    y <- as.numeric(y)
    ylev <- levels(factor(y))
    ynlev <- length(ylev)
    if (ynlev == 0)
        stop("The outcome has no classes/levels!\n")
    ytab <- table(y)
    if (length(ytab) == 1)
        warning("One class of the outcome has no records and will be ignored\n")

    subsample <- function(index, K) {
        permindex <- sample(x=length(index), replace=FALSE, prob=NULL)
        w <- rep(seq_len(K), length.out=length(index))
        out <- list("index"=index[permindex], "which"=w)
        return(out)
    }

    if (K == n) {
        fold <- list("index"=seq_len(n), "which"=seq_len(n))
    } else if (K == 1) {
        fold <- list("index"=seq_len(n), "which"=rep(1, length.out=n))
    } else {
        index <- seq(along = y)
        indexlist <- vector(mode="list", length=ynlev)
        for (l in 1:ynlev) {
            indexlist[[l]] <- index[y == ylev[l]]
        }
        foldlist <- lapply(X=indexlist, FUN=subsample, K=K)
        fold <- list("index"=numeric(0), "which"=numeric(0))
        for (l in 1:ynlev) {
            fold$index <- c(fold$index, foldlist[[l]]$index)
            fold$which <- c(fold$which, foldlist[[l]]$which)
        }
    }
    permkey <- pmatch(x=1:n, table=fold$index)
    ord <- numeric(0)
    for (k in 1:K) {
        ord <- c(ord, fold$index[(fold$which == k)])
    }
    foldkey <- pmatch(x=1:n, table=ord)
    folds <- list(n=n, K=K, perm=fold$index, permkey=permkey, which=fold$which, foldkey=foldkey, seed=seed)

    return(folds)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    endpoints (ind, timemat, probmat, timeval, probval)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

endpoints <- function(ind, timemat, probmat, timeval, probval) {

  N <- nrow(timemat) #N <- nrow(probmat)
  L <- N-length(ind)
  if (!(is.empty(ind))) {
    timemat <- timemat[-ind, , drop=FALSE]
    probmat <- probmat[-ind, , drop=FALSE]
  }
  min.prob.bar <- apply(probmat, 1, min, na.rm=TRUE)
  max.time.bar <- apply(timemat, 1, max, na.rm=TRUE)
  if (is.null(probval) && is.null(timeval)) {
    prob.bar <- rep(NA, L)
    time.bar <- rep(NA, L)
  } else if (!is.null(probval)) {
    prob.bar <- rep(probval, L)
    ind.probmat <- (probmat <= probval)
    ind.probmat[is.na(ind.probmat)] <- FALSE
    time.bar <- numeric(L)
    for (l in 1:L) {
      if (probval >= min.prob.bar[l]) {
        time.bar[l] <- min(timemat[l,which(ind.probmat[l,,drop=TRUE])])
      } else {
        time.bar[l] <- NA
      }
    }
  } else if (!is.null(timeval)) {
    time.bar <- rep(timeval, L)
    ind.timemat <- (timemat >= timeval)
    ind.timemat[is.na(ind.timemat)] <- FALSE
    prob.bar <- numeric(L)
    for (l in 1:L) {
      if (timeval <= max.time.bar[l]) {
        prob.bar[l] <- max(probmat[l,which(ind.timemat[l,,drop=TRUE])])
      } else {
        prob.bar[l] <- NA
      }
    }
  }

  return(list("time.bar"=time.bar,
              "prob.bar"=prob.bar,
              "max.time.bar"=max.time.bar,
              "min.prob.bar"=min.prob.bar))

}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    updatecut (x, fract)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

updatecut <- function(x, fract) {

  p <- dim(x)[2]
  cutpts <- apply(x, 2, "quantile", type=7, probs=fract)

  for (j in 1:p) {
    xunique <- sort(unique(x[, j]))
    if (length(xunique) == 1) {
      cutpts[j] <- min(xunique)
    } else {
      if (isTRUE(all.equal(as.single(cutpts[j]), as.single(min(xunique))))) {
        cutpts[j] <- xunique[2]
      }
    }
  }

  return(cutpts)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    lapply.array (X, rowtrunc=NULL, coltrunc=NULL,
#                                  sub=NULL, fill=NA, MARGIN=1:2, FUN, ...)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

lapply.array <- function (X, rowtrunc=NULL, coltrunc=NULL,
                          sub=NULL, fill=NA, MARGIN=1:2, FUN, ...) {
  x <- list2array(list=X, rowtrunc=rowtrunc, coltrunc=coltrunc, sub=sub, fill=fill)
  return(apply(X=x, MARGIN=MARGIN, FUN=FUN, ...))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    lapply.mat (X, coltrunc=NULL,
#                                sub=NULL, fill=NA, MARGIN=2, FUN, ...)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

lapply.mat <- function (X, coltrunc=NULL,
                        sub=NULL, fill=NA, MARGIN=2, FUN, ...) {
  x <- list2mat(list=X, coltrunc=coltrunc, sub=sub, fill=fill)
  return(apply(X=x, MARGIN=MARGIN, FUN=FUN, ...))
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    list2array (list, rowtrunc=NULL, coltrunc=NULL, sub=NULL, fill=NA)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

list2array <- function (list, rowtrunc=NULL, coltrunc=NULL, sub=NULL, fill=NA) {
  if (!is.empty(list)) {
    if (is.null(sub)) {
      my.list <- list
    } else {
      L <- length(list)
      my.list <- vector(mode="list", length=L)
      for (i in 1:L) {
        my.list[[i]] <- list[[i]][[sub]]
      }
    }
    min.row <- min(sapply(my.list, nrow))
    max.row <- max(sapply(my.list, nrow))
    min.col <- min(sapply(my.list, ncol))
    max.col <- max(sapply(my.list, ncol))
    if (!is.null(coltrunc)) {
      if (coltrunc == "min") {
        adjusted.list <- lapply(my.list, function(x) {x[,1:min.col,drop=FALSE]})
      } else if (coltrunc == "max") {
        adjusted.list <- lapply(my.list, function(x) {cbind(x, matrix(data=fill, nrow=nrow(x), ncol=max.col - ncol(x)))})
      } else {
        adjusted.list <- lapply(my.list, function(x, coltrunc) {if (coltrunc <= ncol(x)) {
                                                                  x[,1:coltrunc,drop=FALSE]
                                                                } else if (coltrunc > ncol(x)) {
                                                                  cbind(x, matrix(data=fill, nrow=nrow(x), ncol=coltrunc - ncol(x)))
                                                                }
                                                                }, coltrunc)
      }
    } else {
        adjusted.list <- lapply(my.list, function(x) {cbind(x, matrix(data=fill, nrow=nrow(x), ncol=max.col - ncol(x)))})
    }
    if (!is.null(rowtrunc)) {
      if (rowtrunc == "min") {
        adjusted.list <- lapply(adjusted.list, function(x) {x[1:min.row,,drop=FALSE]})
      } else if (rowtrunc == "max") {
        adjusted.list <- lapply(adjusted.list, function(x) {rbind(x, matrix(data=fill, nrow=max.row - nrow(x), ncol=ncol(x)))})
      } else {

        adjusted.list <- lapply(my.list, function(x, rowtrunc) {if (rowtrunc <= nrow(x)) {
                                                                  x[1:rowtrunc,,drop=FALSE]
                                                                } else if (rowtrunc > nrow(x)) {
                                                                  rbind(x, matrix(data=fill, nrow=rowtrunc - nrow(x), ncol=ncol(x)))
                                                                }
                                                                }, rowtrunc)
      }
    } else {
        adjusted.list <- lapply(adjusted.list, function(x) {rbind(x, matrix(data=fill, nrow=max.row - nrow(x), ncol=ncol(x)))})
    }
    my.array <- array(data=fill, dim=c(nrow(adjusted.list[[1]]), ncol(adjusted.list[[1]]), length(adjusted.list)))
    for(i in 1:length(adjusted.list)) {
      my.array[,,i] <- adjusted.list[[i]]
    }
  } else {
    my.array <- array(data=fill, dim=c(0,0,0))
  }
  return(my.array)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    list2mat (list, coltrunc=NULL, sub=NULL, fill=NA)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

list2mat <- function (list, coltrunc=NULL, sub=NULL, fill=NA) {
  if (!is.empty(list)) {
    if (is.null(sub)) {
      my.list <- list
    } else {
      L <- length(list)
      my.list <- vector(mode="list", length=L)
      for (i in 1:L) {
        my.list[[i]] <- list[[i]][[sub]]
      }
    }
    min.col <- min(sapply(my.list, length))
    max.col <- max(sapply(my.list, length))
    if (!is.null(coltrunc)) {
      if (coltrunc == "min") {
        adjusted.list <- lapply(my.list, function(x) {x[1:min.col]})
      } else if (coltrunc == "max") {
        adjusted.list <- lapply(my.list, function(x) {c(x, rep(fill, max.col - length(x)))})
      } else {
        adjusted.list <- lapply(my.list, function(x, coltrunc) {if (coltrunc <= length(x)) {
                                                                  x[1:coltrunc]
                                                                } else if (coltrunc > length(x)) {
                                                                  c(x, rep(fill, times=coltrunc - length(x)))
                                                                }
                                                                }, coltrunc)
      }
    } else {
        adjusted.list <- lapply(my.list, function(x) {c(x, rep(fill, max.col - length(x)))})
    }
    my.mat <- do.call(rbind, adjusted.list)
  } else {
    my.mat <- matrix(data=fill, nrow=0, ncol=0)
  }
  return(my.mat)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    cbindlist (list, trunc)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

cbindlist <- function(list, trunc) {
    if (!is.empty(list)) {
        max.row <- max(sapply(list, nrow))
        adjusted.list <- lapply(list, function(x) {rbind(x, matrix(data=NA, nrow=max.row - nrow(x), ncol=ncol(x)))})
        my.mat <- adjusted.list[[1]]
        lcl <- length(adjusted.list)
        if (lcl > 1) {
            for(i in 2:lcl){
                my.mat <- cbind(my.mat, adjusted.list[[i]])
            }
        }
        if (missing(trunc)) {
            trunc <- max.row
        }
        my.mat <- my.mat[1:trunc,,drop=FALSE]
    } else {
        my.mat <- matrix(data=NA, nrow=0, ncol=0)
    }
    return(my.mat)
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    is.empty(x)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

is.empty <- function(x) {
  if (is.vector(x)) {
    if((length(x) == 0) || (x == "")) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if (is.matrix(x) || is.data.frame(x)) {
    return( ((nrow(x) == 0) || (ncol(x) == 0)) )
  } else {
    return( ((length(x) == 0) || (x == "")) )
  }
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    rep.mat (X, margin, times)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

rep.mat <- function (X, margin, times) {
   names(X) <- NULL
   if (missing(margin) || missing(times)) stop("missing arguments")
   if ((class(times) == "numeric") || (class(times) == "integer")){
       times <- as.integer(times)
   } else {
        stop("you must enter an integer or numeric for the repeats")
   }
   if (!((class(X) == "numeric") || (class(X) == "character") || (class(X) == "matrix") || (class(X) == "factor"))) {
        stop("you must enter a matrix, a vector, a scalar, or a factor")
   }
   if (times > 1) {
        if (margin == 2) {
                mat <- X
                for (i in 1:(times-1)) {
                    X <- rbind(X, mat)
                }
                return(X)
        } else if (margin == 1) {
                mat <- X
                for (i in 1:(times-1)) {
                     X <- cbind(X, mat)
                }
                return(X)
        } else {
            stop("you must enter 1 or 2 for the respective direction \"-\" or \"|\"")
        }
   } else if (times == 1) {
        return(X)
   } else {
        return(NULL)
   }
}
###################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    disp (x)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################
disp <- function(x) {
    if (x == "lh") {
        return("LHR")
    } else if (x == "lr") {
        return("LRT")
    } else if (x == "ch") {
        return("CHS")
    } else if (x == "lhr") {
        return("LHR")
    } else if (x == "lrt") {
        return("LRT")
    } else if (x == "cer") {
        return("CER")
    } else if (x == "combined") {
        return("COMBINED")
    } else if (x == "averaged") {
        return("AVERAGED")
    } else if (x == "none") {
        return("NONE")
    }
}
##########################################################################################################################################




##########################################################################################################################################
################
# Usage         :
################
#                    .onAttach (libname, pkgname)
#
################
# Description   :
################
#
################
# Arguments     :
################
#
################
# Values        :
################
#
##########################################################################################################################################

.onAttach <- function(libname, pkgname) {
    SSver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    packageStartupMessage(paste(pkgname, SSver))
    packageStartupMessage("Type PRIMsrc.news() to see new features, changes, and bug fixes\n")
}
##########################################################################################################################################
