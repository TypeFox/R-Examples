fitsurv <- function(parm, Dm, eta) {
  optim(parm, loglik_lamb, gradlik_lamb, method="BFGS", Dm=Dm, eta=eta)$par
}

fitsurv_pw <- function(parm, Dm, eta, breaks) {
  optim(parm, loglik_pw, gradlik_pw, method="BFGS", Dm=Dm, eta=eta, breaks=breaks)$par
}

simoutcome <- function(Xmat, noevent, testtimes, sensitivity, specificity, beta.marker) {
  N <- nrow(Xmat)
  nfeature <- ncol(Xmat)
  blambda <- -log(noevent)/max(testtimes)
  testtimes <- sort(unique(testtimes))
  ntest <- length(testtimes)
  ID <- rep(1:N, each = ntest)
  time <- rep(testtimes, times = N)
  nbiomarker <- length(beta.marker)
  betas <- c(beta.marker, rep(0, nfeature - nbiomarker))
  lambda <- blambda * exp(c(Xmat %*% betas))
  ET <- rexp(N, lambda)
  ET <- ET[ID]
  occur <- time > ET
  probs <- ifelse(occur, sensitivity, 1 - specificity)
  result <- rbinom(length(occur), 1, probs)  
  data <- data.frame(ID, time, result)
  ## NTFP
  afterpos <- function(x) {
    npos <- cumsum(x == 1)
    (npos == 0) | (npos == 1 & x == 1)
  }
  keep <- unlist(tapply(data$result, data$ID, afterpos))
  data <- data[keep, ]
  row.names(data) <- NULL
  data
}

lassofit <- function(Dm, Xmat, minlambda, initsurv, CV.folds=NULL, nlambda = 20, tol=1e-6) {
  J <- ncol(Dm) - 1
  nbeta <- ncol(Xmat)
  parmi <- log(rep(-log(initsurv)/J, J))
  parmi <- c(parmi, rep(0, nbeta))
  maxlam <- maxlambda(Dm, Xmat, parmi[1:J], fitsurv)*1.01
  lambdas <- maxlam*minlambda^((0:nlambda)/nlambda)
  parm <- parmi
  result <- matrix(0, nrow = nbeta+J, ncol=length(lambdas))
  for (i in seq_along(lambdas)) {
    parm <- try(iclasso(Dm, Xmat, parm, lambdas[i], fitsurv, tol))
    if (inherits(parm, "try-error")) break
    result[, i] <- parm
    cat(paste("lambda =", lambdas[i], "finised.\n"))
  }
  if (inherits(parm, "try-error")) {
    cat("RESTART...\n")
    minlambda <- lambdas[i-1]/maxlam
    lambdas <- maxlam*minlambda^((0:nlambda)/nlambda)
    parm <- parmi
    result <- matrix(0, nrow = nbeta+J, ncol=length(lambdas))
    for (i in seq_along(lambdas)) {
      parm <- iclasso(Dm, Xmat, parm, lambdas[i], fitsurv, tol)
      result[, i] <- parm
      cat(paste("lambda =", lambdas[i], "finised.\n"))
    }    
  }
  if (!is.null(CV.folds)) {
    loglik <- function(parm, Dm, Xmat)  -loglikB(parm, Dm, Xmat)
    group <- split(sample(1:nrow(Dm)), 1:CV.folds)
    liks <- matrix(0, nrow = length(group), ncol = length(lambdas))
    for (i in seq_along(group)) {
      id <- group[[i]]
      Dm1 <- Dm[-id, , drop=F]
      Dm2 <- Dm[id, , drop=F]
      Xmat1 <- Xmat[-id, , drop=F]
      Xmat2 <- Xmat[id, , drop=F]
      parm <- parmi
      for (j in seq_along(lambdas)) {
        parm.temp <- try(iclasso(Dm1, Xmat1, parm, lambdas[j], fitsurv, tol))
        if (inherits(parm.temp, "try-error")) {
          liks[i, j] <- NA
        } else {
          parm <- parm.temp
          liks[i, j] <- loglik(parm, Dm2, Xmat2)
        }      
        cat(paste("group", i, "lambdas", j, "finised.\n"))
      }
    }
    liks <- colSums(liks)
  } else {
    liks <- NA
  }
  list(result = result, liks = liks, lambdas = lambdas)
}


lassofit_pw <- function(Dm, Xmat, minlambda, initsurv, breaks, CV.folds=NULL, nlambda = 20, tol=1e-6) {
  J <- ncol(Dm) - 1
  JS <- length(breaks)
  if (breaks[JS] != J-1) stop("invalid breaks")
  nbeta <- ncol(Xmat)
  parmi <- log(rep(-log(initsurv)/J, JS))
  parmi <- c(parmi, rep(0, nbeta))
  maxlam <- maxlambda_pw(Dm, Xmat, parmi[1:JS], breaks, fitsurv_pw)*1.01
  lambdas <- maxlam*minlambda^((0:nlambda)/nlambda)
  parm <- parmi
  result <- matrix(0, nrow = nbeta+JS, ncol=length(lambdas))
  for (i in seq_along(lambdas)) {
    parm <- try(iclasso_pw(Dm, Xmat, parm, breaks, lambdas[i], fitsurv_pw, tol))
    if (inherits(parm, "try-error")) break
    result[, i] <- parm
    cat(paste("lambda =", lambdas[i], "finised.\n"))
  }
  if (inherits(parm, "try-error")) {
    cat("RESTART...\n")
    minlambda <- lambdas[i-1]/maxlam
    lambdas <- maxlam*minlambda^((0:nlambda)/nlambda)
    parm <- parmi
    result <- matrix(0, nrow = nbeta+JS, ncol=length(lambdas))
    for (i in seq_along(lambdas)) {
      parm <- iclasso_pw(Dm, Xmat, parm, breaks, lambdas[i], fitsurv_pw, tol)
      result[, i] <- parm
      cat(paste("lambda =", lambdas[i], "finised.\n"))
    }    
  }
  if (!is.null(CV.folds)) {
    loglik <- function(parm, breaks, Dm, Xmat) {
      par <- parm[1:JS]
      beta <- parm[-(1:JS)]
      eta <- c(Xmat%*%beta)
      -loglik_pw(par, Dm, eta, breaks)
    }
    group <- split(sample(1:nrow(Dm)), 1:CV.folds)
    liks <- matrix(0, nrow = length(group), ncol = length(lambdas))
    for (i in seq_along(group)) {
      id <- group[[i]]
      Dm1 <- Dm[-id, , drop=F]
      Dm2 <- Dm[id, , drop=F]
      Xmat1 <- Xmat[-id, , drop=F]
      Xmat2 <- Xmat[id, , drop=F]
      parm <- parmi
      for (j in seq_along(lambdas)) {
        parm.temp <- try( iclasso_pw(Dm1, Xmat1, parm, breaks, lambdas[j], fitsurv_pw, tol))
        if (inherits(parm.temp, "try-error")) {
          liks[i, j] <- NA
        } else {
          parm <- parm.temp
          liks[i, j] <- loglik(parm, breaks, Dm2, Xmat2)
        }      
        cat(paste("group", i, "lambdas", j, "finised.\n"))
      }
    }
    liks <- colSums(liks)
  } else {
    liks <- NA
  }
  list(result = result, liks = liks, lambdas = lambdas)
}


lassofit_raw <- function(Dm, Xmat, minlambda, initsurv, CV.folds=NULL, nlambda = 20, tol=1e-6) {
  sdv <- Xmat_norm(Xmat)
  J <- ncol(Dm) - 1
  nbeta <- ncol(Xmat)
  parmi <- log(rep(-log(initsurv)/J, J))
  parmi <- c(parmi, rep(0, nbeta))
  maxlam <- maxlambda_raw(Dm, Xmat, sdv, parmi[1:J], fitsurv)*1.01
  lambdas <- maxlam*minlambda^((0:nlambda)/nlambda)
  parm <- parmi
  result <- matrix(0, nrow = nbeta+J, ncol=length(lambdas))
  for (i in seq_along(lambdas)) {
    parm <- try(iclasso_raw(Dm, Xmat, sdv, parm, lambdas[i], fitsurv, tol))
    if (inherits(parm, "try-error")) break
    result[, i] <- parm
    cat(paste("lambda =", lambdas[i], "finised.\n"))
  }
  if (inherits(parm, "try-error")) {
    cat("RESTART...\n")
    minlambda <- lambdas[i-1]/maxlam
    lambdas <- maxlam*minlambda^((0:nlambda)/nlambda)
    parm <- parmi
    result <- matrix(0, nrow = nbeta+J, ncol=length(lambdas))
    for (i in seq_along(lambdas)) {
      parm <- iclasso_raw(Dm, Xmat, sdv, parm, lambdas[i], fitsurv, tol)
      result[, i] <- parm
      cat(paste("lambda =", lambdas[i], "finised.\n"))
    }    
  }
  if (!is.null(CV.folds)) {
    loglik <- function(parm, Dm, Xmat, sdv)  loglik_raw(parm, Dm, Xmat, sdv)
    group <- split(sample(1:nrow(Dm)), 1:CV.folds)
    liks <- matrix(0, nrow = length(group), ncol = length(lambdas))
    for (i in seq_along(group)) {
      id <- group[[i]]
      Dm1 <- Dm[-id, , drop=F]
      Dm2 <- Dm[id, , drop=F]
      Xmat1 <- Xmat[-id, , drop=F]
      Xmat2 <- Xmat[id, , drop=F]
      parm <- parmi
      for (j in seq_along(lambdas)) {
        parm.temp <- try(iclasso_raw(Dm1, Xmat1, sdv, parm, lambdas[j], fitsurv, tol))
        if (inherits(parm.temp, "try-error")) {
          liks[i, j] <- NA
        } else {
          parm <- parm.temp
          liks[i, j] <- loglik(parm, Dm2, Xmat2, sdv)
        }      
        cat(paste("group", i, "lambdas", j, "finised.\n"))
      }
    }
    liks <- colSums(liks)
  } else {
    liks <- NA
  }
  list(result = result, liks = liks, lambdas = lambdas)
}



lassofit_pw_raw <- function(Dm, Xmat, minlambda, initsurv, breaks, CV.folds=NULL, nlambda = 20, tol=1e-6) {
  sdv <- Xmat_norm(Xmat)
  J <- ncol(Dm) - 1
  JS <- length(breaks)
  if (breaks[JS] != J-1) stop("invalid breaks")
  nbeta <- ncol(Xmat)
  parmi <- log(rep(-log(initsurv)/J, JS))
  parmi <- c(parmi, rep(0, nbeta))
  maxlam <- maxlambda_pw_raw(Dm, Xmat, sdv, parmi[1:JS], breaks, fitsurv_pw)*1.01
  lambdas <- maxlam*minlambda^((0:nlambda)/nlambda)
  parm <- parmi
  result <- matrix(0, nrow = nbeta+JS, ncol=length(lambdas))
  for (i in seq_along(lambdas)) {
    parm <- try(iclasso_pw_raw(Dm, Xmat, sdv, parm, breaks, lambdas[i], fitsurv_pw, tol))
    if (inherits(parm, "try-error")) break
    result[, i] <- parm
    cat(paste("lambda =", lambdas[i], "finised.\n"))
  }
  if (inherits(parm, "try-error")) {
    cat("RESTART...\n")
    minlambda <- lambdas[i-1]/maxlam
    lambdas <- maxlam*minlambda^((0:nlambda)/nlambda)
    parm <- parmi
    result <- matrix(0, nrow = nbeta+JS, ncol=length(lambdas))
    for (i in seq_along(lambdas)) {
      parm <- iclasso_pw_raw(Dm, Xmat, sdv, parm, breaks, lambdas[i], fitsurv_pw, tol)
      result[, i] <- parm
      cat(paste("lambda =", lambdas[i], "finised.\n"))
    }    
  }
  if (!is.null(CV.folds)) {
    loglik <- function(parm, breaks, Dm, Xmat, sdv) loglik_pw_raw(parm, breaks, Dm, Xmat, sdv)
    group <- split(sample(1:nrow(Dm)), 1:CV.folds)
    liks <- matrix(0, nrow = length(group), ncol = length(lambdas))
    for (i in seq_along(group)) {
      id <- group[[i]]
      Dm1 <- Dm[-id, , drop=F]
      Dm2 <- Dm[id, , drop=F]
      Xmat1 <- Xmat[-id, , drop=F]
      Xmat2 <- Xmat[id, , drop=F]
      parm <- parmi
      for (j in seq_along(lambdas)) {
        parm.temp <- try(iclasso_pw_raw(Dm1, Xmat1, sdv, parm, breaks, lambdas[j], fitsurv_pw, tol))
        if (inherits(parm.temp, "try-error")) {
          liks[i, j] <- NA
        } else {
          parm <- parm.temp
          liks[i, j] <- loglik(parm, breaks, Dm2, Xmat2, sdv)
        }      
        cat(paste("group", i, "lambdas", j, "finised.\n"))
      }
    }
    liks <- colSums(liks)
  } else {
    liks <- NA
  }
  list(result = result, liks = liks, lambdas = lambdas)
}

bayesfit <- function(Dm, Xmat, b, om1, om2, niter, psample, initsurv, nreport, breaks = NULL, rawX = FALSE) {
  ## standardize Xmat if Xmat is not in raw format
  if (!rawX) {
    means <- matrix(colMeans(Xmat), nrow=nrow(Xmat), ncol=ncol(Xmat), byrow=T)
    sds <- matrix(apply(Xmat, 2, sd), nrow=nrow(Xmat), ncol=ncol(Xmat), byrow=T)
    Xmat <- (Xmat-means)/sds    
  }
  if (is.null(breaks) && !rawX) {
    fit <- bayesmc(Dm, Xmat, b, om1, om2, niter, psample, initsurv, nreport, fitsurv)
  } else if (!is.null(breaks) && !rawX) {
    fit <- bayesmc_pw(Dm, Xmat, breaks, b, om1, om2, niter, psample, initsurv, nreport, fitsurv_pw)
  } else if (is.null(breaks) && rawX) {
    fit <- bayesmc_raw(Dm, Xmat, b, om1, om2, niter, psample, initsurv, nreport, fitsurv)
  } else if (!is.null(breaks) && rawX) {
    fit <- bayesmc_pw_raw(Dm, Xmat, breaks, b, om1, om2, niter, psample, initsurv, nreport, fitsurv_pw)
  }
  fit
}



