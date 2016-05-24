## Modify to use names element of highlod object?
## Fix so call with length(pheno2) == 1 and > 1 give same results; see JoinTestOutputs.
CMSTtests <- function(cross, 
                      pheno1, 
                      pheno2,
                      Q.chr,
                      Q.pos,
                      addcov1 = NULL, 
                      addcov2 = NULL, 
                      intcov1 = NULL, 
                      intcov2 = NULL, 
                      method = c("par", "non.par", "joint", "all"),
                      penalty = c("bic", "aic", "both"),
                      verbose = FALSE)
{
  if (!any(class(cross) == "cross")) 
    stop("Input should have class \"cross\".")

  if(length(pheno2) > 1)
    return(CMSTtestsList(cross, pheno1, pheno2, Q.chr, Q.pos,
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose))
  
  cross.type <- class(cross)[1]

  ## Fit CMST tests
  HkDesignMatrix <- function(qtlo, cross.type) {
    # Create design matrix (Haley-Knott regression)
    nr <- nrow(qtlo$prob[[1]])
    ng <- length(qtlo$prob)
    if(cross.type == "f2"){
      tmp <- 
        unlist(lapply(qtlo$prob, function(x) cbind(x[,1] - x[,3], x[,2])))
      hkm <- matrix(tmp, nr, 2 * ng)
    }
    if(cross.type == "bc"){
      tmp <- unlist(lapply(qtlo$prob, function(x) x[,1] - x[,2]))
      hkm <- matrix(tmp, nr, ng)
    }
    cbind(rep(1, nr), hkm)
  }

  CreateDesignMatrix <- function(cross, chr, pos, addcov.nms, intcov.nms, 
                                 cross.type) {
    # Creates a design matrix
    le.mar <- length(chr)
    n <- nind(cross)
    if(le.mar > 0) {
      qtlo <- makeqtl(cross, chr, pos, what = "prob")
      XG <- HkDesignMatrix(qtlo, cross.type)
    }
    else {
      XG <- matrix(1,n,1)    
    }
    if(!is.null(intcov.nms)) {
      intcov.dat <- data.frame(cross$pheno[, intcov.nms])
      names(intcov.dat) <- intcov.nms
      int.sub.matrix <- model.matrix( as.formula(paste("~", intcov.nms)), 
        intcov.dat)[,-1]
      covs <- data.frame(cross$pheno[, unique(c(addcov.nms, intcov.nms))])
      names(covs) <- unique(c(addcov.nms, intcov.nms))
      form <- as.formula(paste(" ~ ", paste(names(covs), collapse="+")))
      X <- as.matrix(model.matrix(form, data=covs)[,-1])
      if(ncol(XG) > 1) {
        genobyintcov <- XG[,-1] * int.sub.matrix
        X <- cbind(XG, X, genobyintcov)
      }
      else {
        X <- cbind(XG,X)
      }
    }
    else {
      if(!is.null(addcov.nms)){
        covs <- data.frame(cross$pheno[, addcov.nms])
        names(covs) <- addcov.nms
        form <- as.formula(paste(" ~ ", paste(names(covs), collapse = "+")))
        X <- model.matrix(form, data = covs)[,-1]
        if(is.null(dim(X))) {
          X <- as.matrix(X)
        }
        X <- cbind(XG, X)
      }
      else {
        X <- XG
      }
    }
    X
  }

  ParametricIUCMST <- function(Z) {
    pv <- matrix(NA, 4, 4)
    for (i in 1 : 3) {
      for(j in (i + 1) : 4) {
        pv[i, j] <- pnorm(Z[i, j], lower.tail = FALSE)
      }
    }
    pval.1 <- max(pv[1, 2], pv[1, 3], pv[1, 4])
    pval.2 <- max(1 - pv[1, 2], pv[2, 3], pv[2, 4])
    pval.3 <- max(1 - pv[1, 3], 1 - pv[2, 3], pv[3, 4])
    pval.4 <- max(1 - pv[1, 4], 1 - pv[2, 4], 1 - pv[3, 4])
    c(pval.1, pval.2, pval.3, pval.4)
  }

  NonparametricIUCMST <- function(penalty, n, k, vec.logLik) {
    # Computes the non-parametric CMST
    # k: vector of length 4 with the model dimensions
    # vec.logLik: matrix of loglik scores (n by 4)
    kM <- matrix(rep(k, each = n), n, 4)
    if (penalty == "bic") {
      vec.penal <- vec.logLik - 0.5 * kM * log(n)/n
    }
    if (penalty == "aic") {
      vec.penal <- vec.logLik - kM/n
    }
    pv <- matrix(NA, 4, 4)
    for (i in 1 : 3) {
      for (j in (i + 1) : 4) {
        vec.ratio <- vec.penal[, i] - vec.penal[, j]
        counts <- sum(vec.ratio > 0)
        nn <- n - sum(vec.ratio == 0)
        pv[i, j] <- pbinom(counts - 1, nn, 0.5, lower.tail = FALSE)
        pv[j, i] <- pbinom(counts, nn, 0.5, lower.tail = TRUE)
      }
    }
    pval.1 <- max(pv[1, 2], pv[1, 3], pv[1, 4])
    pval.2 <- max(pv[2, 1], pv[2, 3], pv[2, 4])
    pval.3 <- max(pv[3, 1], pv[3, 2], pv[3, 4])
    pval.4 <- max(pv[4, 1], pv[4, 2], pv[4, 3])
    c(pval.1, pval.2, pval.3, pval.4)
  }

  ParametricJointCMST <- function(Z, Cor.hat) {
    z <- min(Z[1, 2], Z[1, 3], Z[1, 4])
    pval.1 <- 1 - pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[1]])
    z <- min(- Z[1, 2], Z[2, 3], Z[2, 4])
    pval.2 <- 1 - pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[2]])
    z <- min(-Z[1, 3], -Z[2, 3], Z[3, 4])
    pval.3 <- 1 - pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[3]])
    z <- min(-Z[1, 4], -Z[2, 4], -Z[3, 4])
    pval.4 <- 1 - pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[4]])
    c(pval.1, pval.2, pval.3, pval.4)
  }

  GetLogLik <- function(cross, y, n, chr, pos, addcov.nms, intcov.nms, 
                        cross.type) {
    X <- CreateDesignMatrix(cross, chr, pos, addcov.nms, intcov.nms, 
                            cross.type)
    dX <- ncol(X)
    qrX <- qr(X)
    b <- qr.coef(qrX, y)
    RSS <- crossprod(y - X %*% b, y - X %*% b)
    log.lik <- as.vector(- (n/2) - (n/2) * log(2 * pi) - (n/2) * log(RSS/n))
    ss <- RSS/n
    vec.log.lik <- dnorm(y, X %*% b, sqrt(ss), log = TRUE)
    list(log.lik = log.lik, vec.log.lik = vec.log.lik, d = dX, RSS = RSS)
  }

  method <- match.arg(method)
  penalty <- match.arg(penalty)

  to.drop <- DropMissing(cross, c(pheno1, pheno2, addcov1, addcov2, 
                         intcov1, intcov2))
  
  if (!is.null(to.drop)) {
    cross <- subset(cross, ind = -to.drop)
  }
  n <- nind(cross)
  y1 <- cross$pheno[, pheno1]
  y2 <- cross$pheno[, pheno2]

  # log.lik.1 #
  tmp.1 <- GetLogLik(cross, y1, n, Q.chr, Q.pos, addcov1, intcov1, 
                     cross.type)
  tmp.2 <- GetLogLik(cross, y2, n, Q.chr, Q.pos, addcov2, intcov2, 
                     cross.type)
  tmp.1g2 <- GetLogLik(cross, y1, n, NULL, NULL, c(addcov1, pheno2), 
                       intcov1, cross.type)
  tmp.2g1 <- GetLogLik(cross, y2, n, NULL, NULL, c(addcov2, pheno1), 
                       intcov2, cross.type)
  tmp.2g1.j <- GetLogLik(cross, y2, n, Q.chr, Q.pos, c(addcov2, pheno1), 
                         intcov2, cross.type)

  TSS1 <- sum((y1 - mean(y1))^2)
  TSS2 <- sum((y2 - mean(y2))^2)
  R2 <- c(1 - (tmp.1$RSS/TSS1), 1 - (tmp.2$RSS/TSS2))

  loglik <- c(tmp.1$log.lik + tmp.2g1$log.lik,
              tmp.2$log.lik + tmp.1g2$log.lik,
              tmp.2$log.lik + tmp.1$log.lik,
              tmp.1$log.lik + tmp.2g1.j$log.lik)  

  model.dim <- 2 + c(tmp.1$d + tmp.2g1$d,
                     tmp.2$d + tmp.1g2$d,
                     tmp.1$d + tmp.2$d,
                     tmp.1$d + tmp.2g1.j$d)

  BICs <- -2 * loglik + model.dim * log(n)
  AICs <- -2 * loglik + 2 * model.dim

  vec.logLik <- cbind(tmp.1$vec.log.lik + tmp.2g1$vec.log.lik,
                      tmp.2$vec.log.lik + tmp.1g2$vec.log.lik,
                      tmp.2$vec.log.lik + tmp.1$vec.log.lik,
                      tmp.1$vec.log.lik + tmp.2g1.j$vec.log.lik)

  vec.LR <- matrix(NA, n, 6, dimnames = list(NULL, 
                   c("12", "13", "14", "23", "24", "34")))

  vec.LR <- matrix(NA, n, 6)
  # i = 1, 12
  # i = 2, 13
  # i = 3, 14
  # i = 4, 23
  # i = 5, 24
  # i = 6, 34 
  ip <- matrix(c(1, 1, 1, 2, 2, 3, 2, 3, 4, 3, 4, 4), 6, 2)
  for (i in 1 : 6) {
    vec.LR[, i] <- vec.logLik[, ip[i, 1]] - vec.logLik[, ip[i, 2]]
  }
  S.hat <- (1 - 1/n) * cov(vec.LR)
  dimnames(S.hat) <- list(dimnames(vec.LR)[[2]], dimnames(vec.LR)[[2]])

  if (method == "all") {
    Sig.hat <- vector(mode = "list", length = 4) 
    signs.2 <- matrix(c(1, -1, -1, -1, 1, 1, -1, 1, 1), 3, 3)
    signs.3 <- matrix(c(1, 1, -1, 1, 1, -1, -1, -1, 1), 3, 3)
    # 1: 12.13.14
    Sig.hat[[1]] <- S.hat[c(1, 2, 3), c(1, 2, 3)]
    # 2: 21.23.24
    Sig.hat[[2]] <- S.hat[c(1, 4, 5), c(1, 4, 5)] * signs.2
    # 3: 31.32.34
    Sig.hat[[3]] <- S.hat[c(2, 4, 6), c(2, 4, 6)] * signs.3
    # 4: 41.42.43
    Sig.hat[[4]] <- S.hat[c(3, 5, 6), c(3, 5, 6)]
    Cor.hat <- vector(mode = "list", length = 4) 
    for (i in 1 : 4) {
      tmp <- is.positive.definite(Sig.hat[[i]])
      if (!tmp) {
        Sig.hat[[i]] <- make.positive.definite(Sig.hat[[i]])
      }
      Cor.hat[[i]] <- cov2cor(Sig.hat[[i]])
      attr(Sig.hat[[i]], "is.positive.definite") <- tmp
    }
    if (penalty == "both") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.BIC <- ParametricIUCMST(Z.bic)
      pvals.np.BIC <- NonparametricIUCMST("bic", n, model.dim, vec.logLik)
      pvals.j.BIC <- ParametricJointCMST(Z.bic, Cor.hat)
      pvals.p.AIC <- ParametricIUCMST(Z.aic)
      pvals.np.AIC <- NonparametricIUCMST("aic", n, model.dim, vec.logLik)
      pvals.j.AIC <- ParametricJointCMST(Z.aic, Cor.hat)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.p.BIC = pvals.p.BIC,
                  pvals.np.BIC = pvals.np.BIC,
                  pvals.j.BIC = pvals.j.BIC,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.p.AIC = pvals.p.AIC,
                  pvals.np.AIC = pvals.np.AIC,
                  pvals.j.AIC = pvals.j.AIC)   
    }
    else if (penalty == "bic") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.BIC <- ParametricIUCMST(Z.bic)
      pvals.np.BIC <- NonparametricIUCMST("bic", n, model.dim, vec.logLik)
      pvals.j.BIC <- ParametricJointCMST(Z.bic, Cor.hat)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.p.BIC = pvals.p.BIC,
                  pvals.np.BIC = pvals.np.BIC,
                  pvals.j.BIC = pvals.j.BIC)         
    }
    else if (penalty == "aic") {
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.AIC <- ParametricIUCMST(Z.aic)
      pvals.np.AIC <- NonparametricIUCMST("aic", n, model.dim, vec.logLik)
      pvals.j.AIC <- ParametricJointCMST(Z.aic, Cor.hat) 
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.p.AIC = pvals.p.AIC,
                  pvals.np.AIC = pvals.np.AIC,
                  pvals.j.AIC = pvals.j.AIC)     
    }
  }
  else if (method == "joint") {
    Sig.hat <- vector(mode = "list", length = 4) 
    signs.2 <- matrix(c(1, -1, -1, -1, 1, 1, -1, 1, 1), 3, 3)
    signs.3 <- matrix(c(1, 1, -1, 1, 1, -1, -1, -1, 1), 3, 3)
    # 1: 12.13.14
    Sig.hat[[1]] <- S.hat[c(1, 2, 3), c(1, 2, 3)]
    # 2: 21.23.24
    Sig.hat[[2]] <- S.hat[c(1, 4, 5), c(1, 4, 5)] * signs.2
    # 3: 31.32.34
    Sig.hat[[3]] <- S.hat[c(2, 4, 6), c(2, 4, 6)] * signs.3
    # 4: 41.42.43
    Sig.hat[[4]] <- S.hat[c(3, 5, 6), c(3, 5, 6)]
    Cor.hat <- vector(mode = "list", length = 4) 
    for (i in 1 : 4) {
      tmp <- is.positive.definite(Sig.hat[[i]])
      if (!tmp) {
        Sig.hat[[i]] <- make.positive.definite(Sig.hat[[i]])
      }
      Cor.hat[[i]] <- cov2cor(Sig.hat[[i]])
      attr(Sig.hat[[i]], "is.positive.definite") <- tmp
    }
    if (penalty == "both") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.j.BIC <- ParametricJointCMST(Z.bic, Cor.hat)
      pvals.j.AIC <- ParametricJointCMST(Z.aic, Cor.hat)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.j.BIC = pvals.j.BIC,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.j.AIC = pvals.j.AIC)   
    }
    else if (penalty == "bic") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.j.BIC <- ParametricJointCMST(Z.bic, Cor.hat)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.j.BIC = pvals.j.BIC)         
    }
    else if (penalty == "aic") {
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.j.AIC <- ParametricJointCMST(Z.aic, Cor.hat) 
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.j.AIC = pvals.j.AIC)     
    }
  }
  else if (method == "par") {
    if (penalty == "both") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.BIC <- ParametricIUCMST(Z.bic)
      pvals.p.AIC <- ParametricIUCMST(Z.aic)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.p.BIC = pvals.p.BIC,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.p.AIC = pvals.p.AIC)
    }
    else if (penalty == "bic") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.BIC <- ParametricIUCMST(Z.bic)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.p.BIC = pvals.p.BIC)     
    }
    else if (penalty == "aic") {
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.p.AIC <- ParametricIUCMST(Z.aic)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.p.AIC = pvals.p.AIC)
    }
  }
  if (method == "non.par") {
    if (penalty == "both") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.np.BIC <- NonparametricIUCMST("bic", n, model.dim, vec.logLik)
      pvals.np.AIC <- NonparametricIUCMST("aic", n, model.dim, vec.logLik)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.np.BIC = pvals.np.BIC,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.np.AIC = pvals.np.AIC)   
    }
    else if (penalty == "bic") {
      Z.bic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (BICs[i] - BICs[j])
          Z.bic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.np.BIC <- NonparametricIUCMST("bic", n, model.dim, vec.logLik)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  BICs = BICs,
                  Z.bic = Z.bic,
                  pvals.np.BIC = pvals.np.BIC)         
    }
    else if (penalty == "aic") {
      Z.aic <- matrix(NA, 4, 4)
      ii <- 1
      for (i in 1 : 3) {
        for (j in (i + 1) : 4) {
          LRt <- -0.5 * (AICs[i] - AICs[j])
          Z.aic[i, j] <- LRt/sqrt(S.hat[ii, ii] * n)
          ii <- ii + 1
        }
      }
      pvals.np.AIC <- NonparametricIUCMST("aic", n, model.dim, vec.logLik)
      out <- list(pheno1 = pheno1,
                  pheno2 = pheno2,
                  n.ind = n,
                  loglik = loglik,  
                  model.dim = model.dim, 
                  R2 = R2,
                  S.hat = S.hat,
                  AICs = AICs,
                  Z.aic = Z.aic,
                  pvals.np.AIC = pvals.np.AIC)     
    }
  }
  out
}
##############################################################################
CMSTtestsList <- function(cross, 
                          pheno1, 
                          pheno2,
                          Q.chr,
                          Q.pos,
                          addcov1 = NULL, 
                          addcov2 = NULL, 
                          intcov1 = NULL, 
                          intcov2 = NULL, 
                          method = c("par", "non.par", "joint", "all"),
                          penalty = c("bic", "aic", "both"),
                          verbose = TRUE)
{
  if(length(pheno2) == 1)
    return(CMSTtestsList(cross, pheno1, pheno2, Q.chr, Q.pos,
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose))

  cross.type <- class(cross)[1]

  ntests <- length(pheno2)
  nms <- paste(pheno1, pheno2, sep = "_")
  pval.nms <- c("pval.1", "pval.2", "pval.3", "pval.4")
  if (penalty == "both") {
    AIC.nms <- c("AIC.1", "AIC.2", "AIC.3", "AIC.4", "z.12", "z.13", "z.14",
                 "z.23", "z.24", "z.34")
    BIC.nms <- c("BIC.1", "BIC.2", "BIC.3", "BIC.4", "z.12", "z.13", "z.14",
                 "z.23", "z.24", "z.34")
    if (method == "all") {
      out <- vector(mode = "list", length = 9)
      names(out) <- c("R2s", "AIC.stats", "BIC.stats", 
                      "pvals.j.BIC", "pvals.p.BIC", "pvals.np.BIC",
                      "pvals.j.AIC", "pvals.p.AIC", "pvals.np.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      for (i in 4 : 9) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)]) 
        out[[3]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[4]][k,] <- aux$pvals.j.BIC
        out[[5]][k,] <- aux$pvals.p.BIC
        out[[6]][k,] <- aux$pvals.np.BIC
        out[[7]][k,] <- aux$pvals.j.AIC
        out[[8]][k,] <- aux$pvals.p.AIC
        out[[9]][k,] <- aux$pvals.np.AIC
        if(verbose)
          cat("pheno2 = ", k, "\n")   
      }
    }
    else if (method == "par") {
      out <- vector(mode = "list", length = 5)
      names(out) <- c("R2s", "AIC.stats", "BIC.stats", "pvals.p.BIC", 
                      "pvals.p.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      for (i in 4 : 5) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)]) 
        out[[3]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[4]][k,] <- aux$pvals.p.BIC
        out[[5]][k,] <- aux$pvals.p.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "non.par") {
      out <- vector(mode = "list", length = 5)
      names(out) <- c("R2s", "AIC.stats", "BIC.stats", "pvals.np.BIC",
                      "pvals.np.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      for (i in 4 : 5) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)]) 
        out[[3]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[4]][k,] <- aux$pvals.np.BIC
        out[[5]][k,] <- aux$pvals.np.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "joint") {
      out <- vector(mode = "list", length = 5)
      names(out) <- c("R2s", "AIC.stats", "BIC.stats", "pvals.j.BIC",
                      "pvals.j.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      for (i in 4 : 5) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)]) 
        out[[3]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[4]][k,] <- aux$pvals.j.BIC
        out[[5]][k,] <- aux$pvals.j.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
  }
  else if (penalty == "bic") {
    BIC.nms <- c("BIC.1", "BIC.2", "BIC.3", "BIC.4", "z.12", "z.13", "z.14",
                 "z.23", "z.24", "z.34")
    if (method == "all") {
      out <- vector(mode = "list", length = 5)
      names(out) <- c("R2s", "BIC.stats", "pvals.j.BIC", "pvals.p.BIC",
                      "pvals.np.BIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      for (i in 3 : 5) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[3]][k,] <- aux$pvals.j.BIC
        out[[4]][k,] <- aux$pvals.p.BIC
        out[[5]][k,] <- aux$pvals.np.BIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "joint") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "BIC.stats", "pvals.j.BIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[3]][k,] <- aux$pvals.j.BIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "par") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "BIC.stats", "pvals.p.BIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[3]][k,] <- aux$pvals.p.BIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "non.par") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "BIC.stats", "pvals.np.BIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, BIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        zb <- t(aux$Z.bic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$BICs, zb[!is.na(zb)])
        out[[3]][k,] <- aux$pvals.np.BIC
        cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
  } 
  else if (penalty == "aic") {
    AIC.nms <- c("AIC.1", "AIC.2", "AIC.3", "AIC.4", "z.12", "z.13", "z.14",
                 "z.23", "z.24", "z.34")
    if (method == "all") {
      out <- vector(mode = "list", length = 5)
      names(out) <- c("R2s", "AIC.stats", "pvals.j.AIC", "pvals.p.AIC",
                      "pvals.np.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      for (i in 3 : 5) {
        out[[i]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      }
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)])
        out[[3]][k,] <- aux$pvals.j.AIC
        out[[4]][k,] <- aux$pvals.p.AIC
        out[[5]][k,] <- aux$pvals.np.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "joint") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "AIC.stats", "pvals.j.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)])
        out[[3]][k,] <- aux$pvals.j.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "par") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "AIC.stats", "pvals.p.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)])
        out[[3]][k,] <- aux$pvals.p.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
    else if (method == "non.par") {
      out <- vector(mode = "list", length = 3)
      names(out) <- c("R2s", "AIC.stats", "pvals.np.AIC")
      out[[1]] <- matrix(NA, ntests, 2, 
                         dimnames = list(nms, c("R2.Y1 ~ Q", "R2.Y2 ~ Q")))
      out[[2]] <- matrix(NA, ntests, 10, dimnames = list(nms, AIC.nms))
      out[[3]] <- matrix(NA, ntests, 4, dimnames = list(nms, pval.nms))
      for(k in 1 : ntests) {
        aux <- CMSTtests(cross, pheno1, pheno2[k], Q.chr, Q.pos, 
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose)
        za <- t(aux$Z.aic)
        out[[1]][k,] <- aux$R2
        out[[2]][k,] <- c(aux$AICs, za[!is.na(za)])
        out[[3]][k,] <- aux$pvals.np.AIC
        if(verbose)
          cat("pheno2 = ", pheno2[k], "\n")   
      }
    }
  } 
  out
}
#########################################################################################
FitAllTests <- function(cross, pheno1, pheno2, Q.chr, Q.pos, verbose = TRUE)
{
  out <- CMSTtests(cross, pheno1, pheno2, Q.chr, Q.pos, 
                     NULL, NULL, NULL, NULL, "all", "both", verbose)

  nms <- pheno2
  ntests <- length(pheno2)
  out$pvals.cit <- matrix(NA, ntests, 2, dimnames = list(nms, c("pval.1", "pval.2")))

  for(k in 1 : ntests) {
    cit.mar <- find.marker(cross, Q.chr, Q.pos)
    LL <- pull.geno(cross)[, cit.mar]
    GG <- cross$pheno[, pheno1]
    TT <- cross$pheno[, pheno2[k]]
    aux2 <- try(CitTests(LL, GG, TT), silent = TRUE)
    if(class(aux2) != "try-error") {
      out$pvals.cit[k,] <- aux2
    }
    if(verbose)
      cat("CIT pheno2 = ", pheno2[k], "\n")   
  }
  out
}
#########################################################################
CitTests <- function(LL, GG, TT)
{
  no.bootstrap <- 50
  ### remove missing values ###
  sel <- (!is.na(LL)) & (!is.na(GG)) & (!is.na(TT))
  dat_f <- as.data.frame(cbind(LL, GG, TT), stringsAsFactors = FALSE)
  dat_f <- dat_f[sel,]
  names(dat_f) <- c("L", "G", "T")
  Lf <- as.factor(dat_f$L)
  dat_f$L <- as.integer(Lf) - 1
  llevels <- as.integer(levels(as.factor(dat_f$L)))
  dfL <- length(llevels) - 1
  pvec <- rep(NA, 4)

  if(dfL == 2){
    dat_f$L1 <- ifelse(dat_f$L == 1,1,0)
    dat_f$L2 <- ifelse(dat_f$L == 2,1,0)
    fit0 <- lm(T ~ 1, data = dat_f)
    fit1 <- lm(T ~ L1 + L2, data = dat_f)
    fit2 <- lm(G ~ T, data = dat_f)
    fit3 <- lm(T ~ G, data = dat_f)
    fit4 <- lm(G ~ T + L1 + L2, data = dat_f)
    fit5 <- lm(T ~ G + L1 + L2, data = dat_f)
    pvec[1] <- anova(fit0, fit1)$"Pr(>F)"[2]
    pvec[2] <- anova(fit2, fit4)$"Pr(>F)"[2]
    pvec[3] <- anova(fit1, fit5)$"Pr(>F)"[2]
    f_ <- anova(fit3, fit5)$F[2]
    fit1G <- lm(G ~ L1 + L2, data = dat_f)
    alg <- summary(fit1G)$coefficients["(Intercept)", 1]
    blg1 <- summary(fit1G)$coefficients["L1", 1]
    blg2 <- summary(fit1G)$coefficients["L2", 1]
    alt <- summary(fit1)$coefficients["(Intercept)", 1]
    blt1 <- summary(fit1)$coefficients["L1", 1]
    blt2 <- summary(fit1)$coefficients["L2", 1]
    dat_f$eG <- resid(fit1G)
    dat_f$eT <- resid(fit1)
    ss <- dim(dat_f)[1]
    fvecr <- rep(NA, no.bootstrap)
    fvecr_r <- rep(NA, no.bootstrap)
    for (i in 1 : no.bootstrap) {
      nni <- trunc(1 + ss*runif(ss, 0, 1)) ;
      dat_f$G_ <- alg + blg1 * dat_f$L1 + blg2 * dat_f$L2 + dat_f$eG[nni]
      fit_0 <- lm(T ~ G_, data = dat_f)
      fit_1 <- lm(T ~ G_ + L1 + L2, data = dat_f)
      fvecr[i] <- anova(fit_0, fit_1)$F[2]
      dat_f$T_ <- alt + blt1 * dat_f$L1 + blt2 * dat_f$L2 + dat_f$eT[nni]
      fit_0 <- lm(G ~ T_, data = dat_f)
      fit_1 <- lm(G ~ T_ + L1 + L2, data = dat_f)
      fvecr_r[i] <- anova(fit_0, fit_1)$F[2]
    }
  }#End dfL == 2

  if(dfL == 1){
    dat_f$L1 <- ifelse(dat_f$L == 1, 1, 0)
    fit0 <- lm(T ~ 1, data = dat_f)
    fit1 <- lm(T ~ L1, data = dat_f)
    fit2 <- lm(G ~ T, data = dat_f)
    fit3 <- lm(T ~ G, data = dat_f)
    fit4 <- lm(G ~ T + L1, data = dat_f)
    fit5 <- lm(T ~ G + L1, data = dat_f)
    pvec[1] <- anova(fit0, fit1)$"Pr(>F)"[2]
    pvec[2] <- anova(fit2, fit4)$"Pr(>F)"[2]
    pvec[3] <- anova(fit1, fit5)$"Pr(>F)"[2]
    f_ <- anova(fit3, fit5)$F[2]
    fit1G <- lm(G ~ L1, data = dat_f)
    alt <- summary(fit1)$coefficients["(Intercept)", 1]
    blt1 <- summary(fit1)$coefficients["L1", 1]
    alg <- summary(fit1G)$coefficients["(Intercept)", 1]
    blg1 <- summary(fit1G)$coefficients["L1", 1]
    dat_f$eG <- resid(fit1G)
    dat_f$eT <- resid(fit1)
    ss <- dim(dat_f)[1]
    fvecr <- rep(NA, no.bootstrap)
    fvecr_r <- rep(NA, no.bootstrap)
    for (i in 1 : no.bootstrap) {
      nni <- trunc(1 + ss*runif(ss, 0, 1)) 
      dat_f$G_ <- alg + blg1 * dat_f$L1 + dat_f$eG[nni]
      fit_0 <- lm(T ~ G_, data = dat_f)
      fit_1 <- lm(T ~ G_ + L1, data = dat_f)
      fvecr[i] <- anova(fit_0, fit_1)$F[2]
      dat_f$T_ <- alt + blt1 * dat_f$L1 + dat_f$eT[nni]
      fit_0 <- lm(G ~ T_, data = dat_f)
      fit_1 <- lm(G ~ T_ + L1, data = dat_f)
      fvecr_r[i] <- anova(fit_0, fit_1)$F[2]
    }
  } #End dfL == 1

  #####F Method
  fvecr <- fvecr[!is.na(fvecr)]
  df1 <- anova(fit3, fit5)$Df[2]
  df2 <- anova(fit3, fit5)$Res.Df[2]
  fncp <- mean(fvecr, na.rm = TRUE) * (df1/df2) * (df2 - df1) - df1
  if(fncp < 0) fncp <- 0
  ######### Transform F to normal
  npvals <- pf(fvecr, df1, df2, ncp = fncp, lower.tail = TRUE)
  nfvecr <- qnorm(npvals)
  npf <- pf(f_, df1, df2, ncp = fncp, lower.tail = TRUE) #Transform observed F
  zf <- qnorm(npf)
  pvec[4] <- pnorm(zf, mean = 0, sd = sd(nfvecr))
  pvalc <- max(pvec)  ###Causal p-value

  #### Reactive p-value
  fit0G <- lm(G ~ 1, data = dat_f)
  pvec1 <- rep(NA, 4)
  pvec1[1] <- anova(fit0G, fit1G)$"Pr(>F)"[2]
  pvec1[2] <- anova(fit3, fit5)$"Pr(>F)"[2]
  pvec1[3] <- anova(fit1G, fit4)$"Pr(>F)"[2]
  f_ <- anova(fit2, fit4)$F[2]
  #####F Method
  fvecr_r <- fvecr_r[!is.na(fvecr_r)]
  df1 <- anova(fit3, fit5)$Df[2]
  df2 <- anova(fit3, fit5)$Res.Df[2]
  fncp <- mean(fvecr_r, na.rm = TRUE) * (df1/df2) * (df2 - df1) - df1
  if(fncp < 0) fncp <- 0
  ######### Transform F to normal
  npvals <- pf(fvecr_r, df1, df2, ncp = fncp, lower.tail = TRUE)
  nfvecr <- qnorm(npvals)
  npf <- pf(f_, df1, df2, ncp = fncp, lower.tail = TRUE) #Transform observed F
  zf <- qnorm(npf)
  pvec1[4] <- pnorm(zf, mean = 0, sd = sd(nfvecr))
  pvalr <- max(pvec1)  ###Reactive p-value
  ###
  c(pvalc, pvalr)
}
##############################################################################
GetCommonQtls <- function(cross, 
                          pheno1, 
                          pheno2,
                          thr = 3,
                          peak.dist = 5,
                          addcov1 = NULL, 
                          addcov2 = NULL, 
                          intcov1 = NULL, 
                          intcov2 = NULL)
{
  if(length(pheno2) > 1 | length(pheno1) > 1)
    stop("pheno1 and pheno2 must have only one trait")
  
  CreateCovMatrix <- function(cross, cov.names) {
    # get covariates data
    if (!is.null(cov.names)) {
      myformula <- formula(paste("~", paste(cov.names, collapse = "+")))
      out <- model.matrix(myformula, cross$pheno)[, -1]
    }
    else {
      out <- NULL
    }
    out
  }
  FindCommonQtls <- function(cross, scanJ, scan1, scan2, thr, peak.dist) {
    # if multiple common QTLs, returns the strongest one from scanJ
    Q <- NA
    ssJ <- summary(scanJ)
    markers1 <- row.names(summary(scan1, thr))
    markers2 <- row.names(summary(scan2, thr))
    chr1 <- scan1[markers1, 1]
    chr2 <- scan2[markers2, 1]
    match.chr <- chr2[match(chr1, chr2, nomatch = 0)]
    if(length(match.chr) > 0) {
      aux1 <- match(match.chr, chr1)
      aux2 <- match(match.chr, chr2)
      markers1 <- markers1[aux1]
      markers2 <- markers2[aux2]
      peak1 <- scan1[markers1, 2]
      peak2 <- scan2[markers2, 2]
      aux3 <- abs(peak1 - peak2)
      aux4 <- which(aux3 <= peak.dist)
      if(length(aux4) > 0){
        cchr <- match.chr[aux4]
        aux5 <- match(cchr, ssJ[, 1])
        aux6 <- which.max(ssJ[aux5, 3])
        Q.chr <- as.numeric(ssJ[aux5[aux6], 1])
        Q.pos <- ssJ[aux5[aux6], 2]
        Q <- find.pseudomarker(cross, Q.chr, Q.pos, "prob")
        Q <- data.frame(Q, Q.chr, Q.pos, stringsAsFactors = FALSE)
      }
    }
    Q
  }

  to.drop <- DropMissing(cross, c(pheno1, pheno2, addcov1, addcov2,
                         intcov1, intcov2))
  if (!is.null(to.drop)) {
    cross <- subset(cross, ind = -to.drop)
  }
  n <- nind(cross)
  y1 <- cross$pheno[, pheno1]
  y2 <- cross$pheno[, pheno2]
  addcov1.M <- CreateCovMatrix(cross, cov.names = addcov1)
  addcov2.M <- CreateCovMatrix(cross, cov.names = addcov2)
  intcov1.M <- CreateCovMatrix(cross, cov.names = intcov1)
  intcov2.M <- CreateCovMatrix(cross, cov.names = intcov2)
  scan1 <- scanone(cross, pheno.col = find.pheno(cross, pheno1), 
                   method = "hk", intcovar = intcov1.M,
                   addcovar = cbind(addcov1.M, intcov1.M))
  scan2 <- scanone(cross, pheno.col = find.pheno(cross, pheno2), 
                   method = "hk", intcovar = intcov2.M, 
                   addcovar = cbind(addcov2.M, intcov2.M))
  scan2g1 <- scanone(cross, pheno.col = find.pheno(cross, pheno2), 
                     method = "hk", intcovar = intcov2.M,
                     addcovar = cbind(addcov2.M, intcov2.M, y1))
  scanJ <- scan1
  scanJ[, 3] <- scan1[, 3] + scan2g1[, 3]
  FindCommonQtls(cross, scanJ, scan1, scan2, thr, peak.dist)
}
##############################################################################
counts <- function(out, alpha, method=c("aic","bic","cit",
  "par.cmst.joint.aic","par.cmst.aic","non.par.cmst.aic",
  "par.cmst.joint.bic","par.cmst.bic","non.par.cmst.bic"))
{
  ###
  get.counts.1 <- function(M, alpha)
  {
    M1 <- sum((M[,1] <= alpha) & (M[,2] > alpha) & (M[,3] > alpha) & (M[,4] > alpha))
    M2 <- sum((M[,1] > alpha) & (M[,2] <= alpha) & (M[,3] > alpha) & (M[,4] > alpha))
    M3 <- sum((M[,1] > alpha) & (M[,2] > alpha) & (M[,3] <= alpha) & (M[,4] > alpha))
    M4 <- sum((M[,1] > alpha) & (M[,2] > alpha) & (M[,3] > alpha) & (M[,4] <= alpha))
    no.call <- nrow(M) - M1 - M2 - M3 - M4
    output <- data.frame(M1,M2,M3,M4,no.call)
  }
  ###
  get.counts.2 <- function(M, alpha)
  {
    M1 <- sum((M[,1] <= alpha) & (M[,2] > alpha))
    M2 <- sum((M[,1] > alpha) & (M[,2] <= alpha))
    M3 <- sum((M[,1] > alpha) & (M[,2] > alpha))
    no.call <- nrow(M) - M1 - M2 - M3
    output <- data.frame(M1,M2,M3,no.call)
  }
  ###
  get.counts.3 <- function(M)
  {
    MM <- t(apply(M,1,rank))
    M1 <- sum(MM[,1]==1)
    M2 <- sum(MM[,2]==1)
    M3 <- sum(MM[,3]==1)
    M4 <- sum(MM[,4]==1)
    no.call <- nrow(MM) - M1 - M2 - M3 - M4
    output <- data.frame(M1,M2,M3,M4,no.call)
  }
  ###
  if(method=="par.cmst.joint.aic")
    output <- get.counts.1(out$pval.par.cmst.joint.AIC, alpha)
  if(method=="par.cmst.aic")
    output <- get.counts.1(out$pval.par.cmst.iu.AIC, alpha)
  if(method=="non.par.cmst.aic")
    output <- get.counts.1(out$pval.non.par.cmst.iu.AIC, alpha)
  if(method=="par.cmst.joint.bic")
    output <- get.counts.1(out$pval.par.cmst.joint.BIC, alpha)
  if(method=="par.cmst.bic")
    output <- get.counts.1(out$pval.par.cmst.iu.BIC, alpha)
  if(method=="non.par.cmst.bic")
    output <- get.counts.1(out$pval.non.par.cmst.iu.BIC, alpha)
  if(method=="cit")
    output <- get.counts.2(out$pval.cit, alpha)
  if(method=="aic")
    output <- get.counts.3(out$AICs)
  if(method=="bic")
    output <- get.counts.3(out$BICs)
  ###
  output
}
##############################################################################
PerformanceSummariesKo <- function(alpha, nms, val.targets, all.orfs,
                                   tests, cis.index)
{
  rnms <- c("aic", "bic", "j.bic", "p.bic", "np.bic", "j.aic", "p.aic", "np.aic", "cit")
  nt <- length(nms)
  TP <- FP <- TN <- FN <- NC <- data.frame(matrix(0, 9, nt))
  tar <- rep(NA, nt)
  names(TP) <- names(FP) <- names(TN) <- names(FN) <- names(NC) <- nms
  row.names(TP) <- row.names(FP) <- row.names(TN) <- row.names(FN) <- row.names(NC) <- rnms
  Causal <- NotCausal <- vector(mode = "list", length = 9)
  for (k in 1 : nt) {
      aux <- which(tests[[11]][,1] == nms[k])
      aux.nms <- tests[[11]][aux, 2]
      tar[k] <- length(aux)
      for (i in 2 : 3) {
        aux.rank <- apply(tests[[i]][aux, 1:4, drop = F], 1, rank)
        aux.best <- apply(aux.rank, 2, function(x) which.min(x))
        aux.index <- as.numeric(which(aux.best == 1))
        Causal[[i - 1]] <- aux.nms[aux.index]
        NotCausal[[i - 1]] <- aux.nms[-aux.index]
      }
      for (i in 4 : 9) {
        Causal[[i - 1]] <- aux.nms[which(tests[[i]][aux, 1] <= alpha)]
        NotCausal[[i - 1]] <-
          aux.nms[c(which(tests[[i]][aux, 2] <= alpha),
                    which(tests[[i]][aux, 3] <= alpha),
                    which(tests[[i]][aux, 4] <= alpha))]
      }
      Causal[[9]] <-
        aux.nms[which(tests[[10]][aux, 1] <= alpha & tests[[10]][aux, 2] > alpha)]
      NotCausal[[9]] <-
        aux.nms[c(which(tests[[10]][aux, 1] > alpha & tests[[10]][aux, 2] <= alpha),
                  which(tests[[10]][aux, 1] >= alpha & tests[[10]][aux, 2] >= alpha))]
      val <- val.targets[[match(nms[k], names(val.targets))]]
      not.val <- all.orfs[-match(unique(c(nms[k], val)), all.orfs)]
      for (i in 1 : 9) {
        TP[i, k] <- length(!is.na(intersect(Causal[[i]], val)))
        FP[i, k] <- length(!is.na(intersect(Causal[[i]], not.val)))
        TN[i, k] <- length(!is.na(intersect(NotCausal[[i]], not.val)))
        FN[i, k] <- length(!is.na(intersect(NotCausal[[i]], val)))
      }
      for (i in 4 : 9) {
        NC[i - 1, k] <- length(c(which(tests[[i]][aux, 1] > alpha),
                             which(tests[[i]][aux, 2] > alpha),
                             which(tests[[i]][aux, 3] > alpha),
                             which(tests[[i]][aux, 4] > alpha)))
      }
      NC[9, k] <- length(c(which(tests[[10]][aux, 1] < alpha),
                           which(tests[[10]][aux, 2] < alpha)))
  }
  tp <- apply(TP, 1, sum)
  fp <- apply(FP, 1, sum)
  tn <- apply(TN, 1, sum)
  fn <- apply(FN, 1, sum)
  nc <- apply(NC, 1, sum)
  prec <- tp/(tp + fp)
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  overall.1 <- data.frame(prec, tp, fp, tpr, fpr, tn, fn, nc)
  tp <- apply(TP[, cis.index], 1, sum)
  fp <- apply(FP[, cis.index], 1, sum)
  tn <- apply(TN[, cis.index], 1, sum)
  fn <- apply(FN[, cis.index], 1, sum)
  nc <- apply(NC[, cis.index], 1, sum)
  prec <- tp/(tp + fp)
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  overall.2 <- data.frame(prec, tp, fp, tpr, fpr, tn, fn, nc)
  list(overall.1, overall.2, tar)
}

##############################################################################
performance.summaries.cmst <- function(out, model, alpha=0.05, method)
{
  ntests <- nrow(out[[1]])
  if((model=="A") || (model=="B")){
    ct <- counts(out, alpha, method)
    tp <- ct[1,1]
    fp <- ct[1,2]+ct[1,3]+ct[1,4]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)
  } 
  if(model=="C"){
    ct <- counts(out, alpha, method)
    tp <- ct[1,4]
    fp <- ct[1,1]+ct[1,2]+ct[1,3]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp) 
  }
  if(model=="D"){
    ct <- counts(out, alpha, method)
    tp <- ct[1,3]
    fp <- ct[1,1]+ct[1,2]+ct[1,4]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)  
  }
  if(model=="E"){
    ct <- counts(out, alpha, method)
    tp <- ct[1,4]
    fp <- ct[1,1]+ct[1,2]+ct[1,3]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)   
  } 
  data.frame(tp, fp, power, type1.err, prec)
}
##############################################################################
performance.summaries.cit <- function(out, model, alpha=0.05)
{
  ntests <- nrow(out[[1]])
  if((model=="A") || (model=="B")){
    ct <- counts(out, alpha, method="cit")
    tp <- ct[1,1]
    fp <- ct[1,2]+ct[1,3]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)
  } 
  if((model=="C") || (model=="D") || (model=="E")){
    ct <- counts(out, alpha, method="cit")
    tp <- ct[1,3]
    fp <- ct[1,1]+ct[1,2]
    power <- tp/ntests
    type1.err <- fp/ntests
    prec <- tp/(tp+fp)   
  }
  data.frame(tp, fp, power, type1.err, prec)
}
##############################################################################
get.power.type1.prec.matrix <- function(out, models, alpha)
{
  n <- length(alpha)
  Power <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
    "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
    as.character(alpha)))
  Type1 <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
    "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
    as.character(alpha)))
  Prec <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
    "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
    as.character(alpha)))
  for(k in 1:n){
    outs <- array(NA, c(9,2,5), dimnames=list(c("aic","par.joint.aic","par.aic",
      "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
      c("TP","FP"),models))
    for(i in 1:5){
      outs[1,1:2,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="aic")[1:2])
      outs[2,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.joint.aic")[1:2])
      outs[3,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.aic")[1:2])
      outs[4,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="non.par.cmst.aic")[1:2])
      outs[5,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="bic")[1:2])
      outs[6,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.joint.bic")[1:2])
      outs[7,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.bic")[1:2])
      outs[8,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="non.par.cmst.bic")[1:2])
      outs[9,,i] <- as.numeric(performance.summaries.cit(out[[i]], model=models[i], alpha[k])[1:2])
    }
    all <- matrix(NA, 9, 2, dimnames=list(c("aic","par.joint.aic","par.aic",
      "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
      c("TP","FP")))
    for(i in 1:9)
      all[i,] <- apply(outs[i,1:2,],1,sum)
    Prec[,k] <- all[,1]/apply(all,1,sum)
    Power[,k] <- all[,1]/5000
    Type1[,k] <- all[,2]/5000
    print(k)
  }
  list(Power=Power, Type1=Type1, Prec=Prec)
}
##############################################################################
## without model C
get.power.type1.prec.matrix.2 <- function(out, models, alpha)
{
  n <- length(alpha)
  Power <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
    "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
    as.character(alpha)))
  Type1 <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
    "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
    as.character(alpha)))
  Prec <- matrix(NA,9,n, dimnames=list(c("aic","par.joint.aic","par.aic",
    "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
    as.character(alpha)))
  for(k in 1:n){
    outs <- array(NA, c(9,2,5), dimnames=list(c("aic","par.joint.aic","par.aic",
      "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
      c("TP","FP"),models))
    for(i in 1:5){
      outs[1,1:2,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="aic")[1:2])
      outs[2,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.joint.aic")[1:2])
      outs[3,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.aic")[1:2])
      outs[4,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="non.par.cmst.aic")[1:2])
      outs[5,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="bic")[1:2])
      outs[6,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.joint.bic")[1:2])
      outs[7,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="par.cmst.bic")[1:2])
      outs[8,,i] <- as.numeric(performance.summaries.cmst(out[[i]], model=models[i], alpha[k], method="non.par.cmst.bic")[1:2])
      outs[9,,i] <- as.numeric(performance.summaries.cit(out[[i]], model=models[i], alpha[k])[1:2])
    }
    all <- matrix(NA, 9, 2, dimnames=list(c("aic","par.joint.aic","par.aic",
      "non.par.aic","bic","par.joint.bic","par.bic","non.par.bic","cit"),
      c("TP","FP")))
    for(i in 1:9)
      all[i,] <- apply(outs[i,1:2,-3],1,sum)
    Prec[,k] <- all[,1]/apply(all,1,sum)
    Power[,k] <- all[,1]/4000
    Type1[,k] <- all[,2]/4000
    print(k)
  }
  list(Power=Power, Type1=Type1, Prec=Prec)
}
##############################################################################
CombineTests <- function(comap, file)
{
  reg.nms <- names(comap)
  out <- NULL
  join.out <- list()
  for (k in 1 : length(comap)) {
    load(paste(file, reg.nms[k], "Rdata", sep="."))
    join.out[[k]] <- out
  }
  names(join.out) <- reg.nms
  join.out
}
JoinTestOutputs <- function(comap, tests, file = NULL)
{
  if(!is.null(file) & missing(tests))
    tests <- CombineTests(comap, file)
  
  reg.nms <- names(comap)
  join.out <- tests[[1]]
  ## Add extra element to join.out: phenos.
  join.out$phenos <- cbind(rep(reg.nms[1], length(comap[[1]])), comap[[1]])

  for (k in 2 : length(comap)) {
    out <- tests[[k]]
    if (length(out) != 10) { ## CMSTtests if length(pheno2) == 1
      tmp <- t(out$Z.aic)
      AIC.stats <- c(out$AICs, tmp[!is.na(tmp)])
      tmp <- t(out$Z.bic)
      BIC.stats <- c(out$BICs, tmp[!is.na(tmp)])      
      out <- list(R2s = out$R2, AIC.stats = AIC.stats, BIC.stats = BIC.stats, 
                  pvals.j.BIC = out$pvals.j.BIC, pvals.p.BIC = out$pvals.p.BIC, 
                  pvals.np.BIC = out$pvals.np.BIC,  pvals.j.AIC = out$pvals.j.AIC, 
                  pvals.p.AIC = out$pvals.p.AIC,  pvals.np.AIC = out$pvals.np.AIC,
                  pvals.cit = out$pvals.cit)
    }

    for (i in 1 : 10) {
      join.out[[i]] <- rbind(join.out[[i]], out[[i]])
    }
    join.out[[11]] <- 
      rbind(join.out[[11]], cbind(rep(reg.nms[k], length(comap[[k]])), comap[[k]]))
  }
 join.out
}
##############################################################################
p.adjust.np <- function(tests, method = "BH")
{
  for(test.name in paste("pvals.np", c("AIC","BIC"), sep = "."))
    for (i in seq(ncol(tests[[test.name]])))
      tests[[test.name]][, i] <- p.adjust(tests[[test.name]][, i], method = method)

  tests

}
##############################################################################
GetCis <- function(x, window = 10) {
  xx <- x[x[, 2] == x[, 4],]
  xx <- xx[abs(xx[, 3] - xx[, 5]) <= window, ]
  index <- match(xx[, 1], x[, 1])
  list(cis.reg = xx, cis.index = index) 
}
##############################################################################
GetCisCandReg <- function(highobj, cand.reg, lod.thr = NULL)
{
  cand.names <- as.character(cand.reg[, 1])
  
  ## Restrict to being on same chromosome. This is fragile.
  cand.reg <- cand.reg[cand.reg[, 2] == cand.reg[, 4],]

  ## Restrict to LOD above lod.thr.
  highobj <- highlod.thr(highobj, lod.thr)
  
  chr.pos <- highobj$chr.pos

  ## Subset highlod to those phenos in cand.reg.
  pheno.cols <- unique(highobj$highlod$phenos)
  m.pheno <- match(as.character(cand.reg[,1]), highobj$names[pheno.cols])
  if(any(is.na(m.pheno)))
    stop("cannot find cand.reg traits in highobj$names")
  m.pheno <- pheno.cols[m.pheno]
  highlod <- highobj$highlod[highobj$highlod$phenos %in% m.pheno, ]
  
  ## Get start and end for each pheno. NB: may include multiple chr.
  h.index <- cumsum(table(highlod$phenos))
  h.index <- cbind(start = 1 + c(0, h.index[-length(h.index)]), end = h.index)
  ## Now get in right order.
  m.pheno <- order(m.pheno)
  h.index[m.pheno,] <- h.index
  
  ## Find lower and upper position around peak.
  tmpfn <- function(x, highlod, chr.pos) {
    h <- highlod[x[1]:x[2],, drop = FALSE]
    ## Only look at chr with peak LOD.
    wh <- which.max(h$lod)
    wh <- range(which(chr.pos$chr[h$row] == chr.pos$chr[h$row[wh]]))
    ## Could have non-contiguous regions. Don't sweat it for now.
    chr.pos$pos[h$row[wh]]
  }
  peak.pos <- t(apply(h.index, 1, tmpfn, highlod, chr.pos))
  dimnames(peak.pos)[[2]] <- c("peak.pos.lower", "peak.pos.upper")

  out <- data.frame(cand.reg, peak.pos)
  is.cis <- (out$phys.pos >= out$peak.pos.lower &
             out$phys.pos <= out$peak.pos.upper)
  ## Keep cis traits, but leave off peak.chr (since it == phys.chr).
  out <- out[is.cis, -4, drop = FALSE]
  if(nrow(out))
    attr(out, "cis.index") <- match(out[, 1], cand.names)
  out
}
##############################################################################
PerformanceSummariesKo <- function(alpha, nms, val.targets, all.orfs, 
                                   tests, cis.index)
{
  ## Unclear what this does. Part of KO data analysis.
  ## tests added.
  
  rnms <- c("aic", "bic", "j.bic", "p.bic", "np.bic", "j.aic", "p.aic", "np.aic", "cit")
  nt <- length(nms)
  TP <- FP <- TN <- FN <- NC <- data.frame(matrix(0, 9, nt))
  tar <- rep(NA, nt)
  names(TP) <- names(FP) <- names(TN) <- names(FN) <- names(NC) <- nms
  row.names(TP) <- row.names(FP) <- row.names(TN) <- row.names(FN) <- row.names(NC) <- rnms
  Causal <- NotCausal <- vector(mode = "list", length = 9)
  length.intersect <- function(a, b) {
    ints <- intersect(a, b)
    if(is.null(ints))
      0
    else
      length(!is.na(ints))
  }
  for (k in 1 : nt) {
      aux <- which(tests[[11]][,1] == nms[k])
      aux.nms <- tests[[11]][aux, 2]
      tar[k] <- length(aux)
      for (i in 2 : 3) {
        aux.rank <- apply(tests[[i]][aux, 1:4, drop = F], 1, rank)
        aux.best <- apply(aux.rank, 2, function(x) which.min(x))
        aux.index <- as.numeric(which(aux.best == 1))
        Causal[[i - 1]] <- aux.nms[aux.index] 
        NotCausal[[i - 1]] <- aux.nms[-aux.index] 
      }  
      for (i in 4 : 9) {
        Causal[[i - 1]] <- aux.nms[which(tests[[i]][aux, 1] <= alpha)]
        NotCausal[[i - 1]] <- 
          aux.nms[c(which(tests[[i]][aux, 2] <= alpha),
                    which(tests[[i]][aux, 3] <= alpha),
                    which(tests[[i]][aux, 4] <= alpha))] 
      }  
      Causal[[9]] <- 
        aux.nms[which(tests[[10]][aux, 1] <= alpha & tests[[10]][aux, 2] > alpha)]
      NotCausal[[9]] <- 
        aux.nms[c(which(tests[[10]][aux, 1] > alpha & tests[[10]][aux, 2] <= alpha),
                  which(tests[[10]][aux, 1] >= alpha & tests[[10]][aux, 2] >= alpha))]
      val <- val.targets[[match(nms[k], names(val.targets))]]
      not.val <- all.orfs[-match(unique(c(as.character(nms[k]), val)), all.orfs)]
      for (i in 1 : 9) {
        TP[i, k] <- length.intersect(Causal[[i]], val)
        FP[i, k] <- length.intersect(Causal[[i]], not.val)
        TN[i, k] <- length.intersect(NotCausal[[i]], not.val)
        FN[i, k] <- length.intersect(NotCausal[[i]], val)
      }
      for (i in 4 : 9) {
        NC[i - 1, k] <- length(c(which(tests[[i]][aux, 1] > alpha),
                             which(tests[[i]][aux, 2] > alpha),
                             which(tests[[i]][aux, 3] > alpha),
                             which(tests[[i]][aux, 4] > alpha)))
      }
      NC[9, k] <- length(c(which(tests[[10]][aux, 1] < alpha),
                           which(tests[[10]][aux, 2] < alpha)))
  }
  tp <- apply(TP, 1, sum)
  fp <- apply(FP, 1, sum)
  tn <- apply(TN, 1, sum)
  fn <- apply(FN, 1, sum)
  nc <- apply(NC, 1, sum)
  prec <- tp/(tp + fp)
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  overall.1 <- data.frame(prec, tp, fp, tpr, fpr, tn, fn, nc)
  tp <- apply(TP[, cis.index, drop = FALSE], 1, sum)
  fp <- apply(FP[, cis.index, drop = FALSE], 1, sum)
  tn <- apply(TN[, cis.index, drop = FALSE], 1, sum)
  fn <- apply(FN[, cis.index, drop = FALSE], 1, sum)
  nc <- apply(NC[, cis.index, drop = FALSE], 1, sum)
  prec <- tp/(tp + fp)
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  overall.2 <- data.frame(prec, tp, fp, tpr, fpr, tn, fn, nc)
  list(overall.1, overall.2, tar)
}
##############################################################################
PrecTpFpMatrix <- function(alpha, val.targets, all.orfs, tests, cand.reg, cis.cand.reg)
{
  nms <- as.character(cand.reg[,1])
  cis.index <- attr(cis.cand.reg, "cis.index")
  
  le <- length(alpha)
  Prec1 <- Tp1 <- Fp1 <- matrix(NA, 9, le, 
    dimnames=list(c("aic", "bic", "j.bic", "p.bic", "np.bic", "j.aic", 
                    "p.aic", "np.aic", "cit"), as.character(alpha)))
  Prec2 <- Tp2 <- Fp2 <- Prec1
  for(i in 1:le){
    aux <- PerformanceSummariesKo(alpha = alpha[i], nms,
                                  val.targets = val.targets, 
                                  all.orfs = all.orfs, 
                                  tests = tests,
                                  cis.index = cis.index)
    Prec1[,i] <- round(aux[[1]][,1], 2) 
    Prec2[,i] <- round(aux[[2]][,1], 2)
    Tp1[,i] <- aux[[1]][, 2]
    Tp2[,i] <- aux[[2]][, 2]
    Fp1[,i] <- aux[[1]][, 3]
    Fp2[,i] <- aux[[2]][, 3]
  }
  list(Prec1 = Prec1,
       Prec2 = Prec2,
       Tp1 = Tp1,
       Tp2 = Tp2,
       Fp1 = Fp1,
       Fp2 = Fp2)
}
##############################################################################
CreateTraitsLodInt <- function(scan, annot, traits, lod.thr, drop = 1.5)
{
  traits <- unique(traits)
  n <- length(traits)
  out <- data.frame(matrix(NA, n, 7))
  names(out) <- c("gene", "chr", "phys.pos", "lower.pos", "peak.pos",
                  "upper.pos", "lod")
  nms <- names(scan)[-c(1, 2)]
  for (i in 1 : n) {
    ii <- match(traits[i], annot[, 1])
    out[i, 1:3] <- annot[ii, c(1, 3, 5)]
    trait.chr <- annot[ii, 3]
    trait.pos <- annot[ii, 5]
    if (!is.na(trait.pos)) {
      peak <- max(scan[scan[, 1] == trait.chr, traits[i]])
      if(peak >= lod.thr){
        trait.index <- match(traits[i], nms)
        sscan <- scan[, c(1, 2, trait.index + 2)]
        lod.interval <- lodint(sscan, chr = trait.chr, drop)
        lb <- lod.interval[1, 2]
        ub <- lod.interval[3, 2]
        out[i, 4] <- lb
        out[i, 5] <- lod.interval[2, 2]
        out[i, 6] <- ub
        out[i, 7] <- peak
      }
    }     
    cat(" ", i, "\n")
  }
  subset(out, !is.na(out[, 4]))
}
##############################################################################
GetCandReg <- function(highobj, annot, traits)
{
  ## currently this only gets max over genome; want max by chr, yes?
  traits <- unique(traits)
  n <- length(traits)
  out <- data.frame(matrix(NA, n, 6))
  names(out) <- c("gene", "phys.chr", "phys.pos", "peak.chr", "peak.pos",
                  "peak.lod")

  ## Get annotation of trait name and physical chromosome and position (in cM).
  m <- match(traits, annot[,1])
  out[, 1:3] <- annot[m, c(1,3,5)]

  ## Get LOD peak information
  m <- !is.na(out[,3])

  pheno.cols <- match(traits[m], highobj$names)
  if(any(is.na(pheno.cols)))
     stop("some traits do not have scans")

  chr.pos <- highobj$chr.pos
  highlod <- highobj$highlod[highobj$highlod$phenos %in% pheno.cols,]
  tmp <- cumsum(table(highlod[,"phenos"]))
  tmp <- c(0, tmp[-length(tmp)])
  peak.index <- tmp + tapply(highlod$lod, highlod$phenos, which.max)

  ## now relate to lod, chr, pos, but get order right with traits
  m <- match(highobj$names[highlod[peak.index, "phenos"]], as.character(out[,1]))
  if(any(is.na(m)))
    stop("cannot match highlod with pheno names")
  
  out[m, 6] <- highlod[peak.index, "lod"]
  out[m, 4] <- chr.pos[highlod[peak.index, "row"], "chr"]
  out[m, 5] <- chr.pos[highlod[peak.index, "row"], "pos"]

  out[!is.na(out[,4]),, drop = FALSE]
}
##############################################################################
GetCoMappingTraits <- function(highobj, cand.reg)
{
  chr.pos <- highobj$chr.pos
  chrs <- levels(chr.pos$chr)
  chr <- ordered(cand.reg$peak.chr, chrs)
  phys <- cbind(phenos = match(as.character(cand.reg[,1]), highobj$names),
                chr = unclass(chr), pos = cand.reg$peak.pos)

  in.range <- function(pos, x.pos) {
    r <- range(pos)
    r[1] <= x.pos & r[2] >= x.pos
  }
  ## Find traits that are co-mapping.
  find.comap <- function(x, highlod, chr.pos, chrs, all.traits) {
    ## Traits must map to same chromosome.
    h <- highlod[chrs[x[2]] == chr.pos$chr[highlod$row],, drop = FALSE]
    ## And peaks must be in range.
    h <- tapply(chr.pos$pos[h$row], h$phenos, in.range, x[3])
    ## But have to remove trait from its list.
    h <- as.numeric(names(h[h]))
    h <- h[-match(x[1], h)]
    all.traits[h]
  }

  ## This list is too restrictive compared with earlier list of Elias.
  ## Try re-running deprecated code using scan.orf to compare.

  out <- apply(phys, 1, find.comap, highobj$highlod, chr.pos, chrs, highobj$names)
  if(is.matrix(out))
    out <- as.data.frame(out)
  names(out) <- as.character(cand.reg[,1])

  lapply(out, as.character)
}
