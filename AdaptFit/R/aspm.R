########## R-function: aspm ##########

# Modified version of spm function to take 
# estimated variance of random effects into account

# Last changed: 30 JUL 2007


"aspm" <-
  function (spm.info, random = NULL, group = NULL, family = "gaussian", 
            spar.method = "REML", omit.missing = NULL, Si.b = NULL, weights = NULL,correlation=NULL,control=NULL) 
{
  require("nlme")
  
  random.info <- NULL
  if (!is.null(random)) 
    random.info <- random.read(random, group)
  spm.info <- c(spm.info, list(random = random.info))
  if (!is.null(unlist(spm.info$pen$spar))) {
    if (any(unlist(spm.info$pen$spar) == 0)) 
      stop("zero smoothing parameters not supported in current version.")
  }
  design.info <- spmDesign(spm.info)
  X <- design.info$X
  Z <- design.info$Z

#if the variance of random effects Si.b is given, spline basis matrix is standardized

  if (!is.null(Si.b)) 
    Z <- t(t(Z) * sqrt(Si.b))
  y <- spm.info$y
  trans.mat <- design.info$trans.mat
  spm.info <- c(spm.info, list(trans.mat = trans.mat))
  block.inds <- design.info$block.inds
  re.block.inds <- design.info$re.block.inds
  col.ones <- rep(1, nrow(X))
  auto.spar.select <- FALSE
  if ((is.null(spm.info$pen)) & (is.null(spm.info$krige))) 
    auto.spar.select <- TRUE
  if (!is.null(spm.info$pen)) {
    auto.spar <- 0
    for (j in 1:length(spm.info$pen$name)) auto.spar <- (auto.spar + 
                                                         (is.null(spm.info$pen$spar[[j]]) & (spm.info$pen$adf[[j]] == 
                                                                                             "miss")))
    if (auto.spar > 0) 
      auto.spar.select <- TRUE
  }
  if (!is.null(spm.info$krige)) {
    if ((is.null(spm.info$krige$spar)) & (spm.info$krige$adf[[1]] == 
                                          "miss")) 
      auto.spar.select <- TRUE
  }
  if (auto.spar.select == FALSE) {
    if (!is.null(spm.info$lin)) 
      num.lin <- ncol(as.matrix(spm.info$lin$x))
    if (is.null(spm.info$lin)) 
      num.lin <- 0
    compon.num <- 1
    if (!is.null(spm.info$pen)) {
      basis.type <- spm.info$pen$basis
      for (j in 1:ncol(as.matrix(spm.info$pen$x))) {
        if (!is.null(spm.info$pen$adf[[j]]) & (spm.info$pen$adf[[j]] != "miss")) {
          deg.val <- spm.info$pen$degree[j]
          stt.ind <- block.inds[[compon.num + num.lin + 1]][1]
          ncol.Xj <- length(block.inds[[compon.num + num.lin + 1]])
          ncol.Xj <- ncol.Xj - length(re.block.inds[[compon.num]])
          Xj <- X[, stt.ind:(stt.ind + ncol.Xj - 1)]
          Xj <- cbind(col.ones, Xj)
          Zj <- Z[, re.block.inds[[compon.num]]]
          adf.val <- spm.info$pen$adf[[j]]
          if (family == "gaussian") 
            spar.val <- df.to.spar(adf.val + 1, Xj, Zj)
          if (!family == "gaussian") 
            spar.val <- glm.df.to.spar(adf.val + 1, y,  Xj, Zj, family)
          if (basis.type == "trunc.poly") 
            spm.info$pen$spar[[j]] <- exp(log(spar.val)/(2 *  deg.val))
          else spm.info$pen$spar[[j]] <- exp(log(spar.val)/deg.val)
          compon.num <- compon.num + 1
        }
      }
    }
    if (!is.null(spm.info$krige)) {
      if ((!is.null(spm.info$krige$adf)) & (spm.info$krige$adf[[1]] != 
                                            "miss")) {
        deg.val <- spm.info$krige$degree
        stt.ind <- block.inds[[compon.num + num.lin + 1]][1]
        ncol.Xj <- length(block.inds[[compon.num + num.lin + 1]])
        ncol.Xj <- ncol.Xj - length(re.block.inds[[compon.num]])
        Xj <- X[, stt.ind:(stt.ind + ncol.Xj - 1)]
        Xj <- cbind(col.ones, Xj)
        Zj <- Z[, re.block.inds[[compon.num]]]
        adf.val <- spm.info$krige$adf[[1]]
        if (family == "gaussian") 
          spar.val <- df.to.spar(adf.val + 1, Xj, Zj)
        if (!family == "gaussian") 
          if (is.null(spm.info$off.set)) 
            spar.val <- glm.df.to.spar(adf.val + 1, y, Xj, Zj, family)
          else spar.val <- glm.df.to.spar(adf.val + 1, y - spm.info$off.set, Xj, Zj, family)
        spm.info$krige$spar <- exp(log(spar.val)/deg.val)
      }
    }
    diag.G <- NULL
    if (!is.null(spm.info$pen)) 
      for (j in 1:ncol(as.matrix(spm.info$pen$x))) {
        deg.val <- spm.info$pen$degree[j]
        spar.val <- spm.info$pen$spar[[j]]
        num <- length(spm.info$pen$knots[[j]])
        if (basis.type == "trunc.poly") 
          diag.G <- c(diag.G, rep(1/(exp((2 * deg.val) *  log(spar.val))), num))
        else diag.G <- c(diag.G, rep(1/(exp((deg.val) *   log(spar.val))), num))
      }
    if (!is.null(spm.info$krige)) {
      spar.val <- spm.info$krige$spar
      num.knots <- nrow(spm.info$krige$knots)
      diag.G <- c(diag.G, rep((1/spar.val^2), num.knots))
    }
  }
  group.vec <- col.ones
  assign("group.vec.Handan", group.vec, pos = 1)
  assign("X.Declan", X, pos = 1)
  if (!is.null(Z)) {
    assign("Z.Jaida", Z, pos = 1)
    data.fr <- groupedData(y ~ -1 + X.Declan | group.vec.Handan,  data = data.frame(y, X.Declan, Z.Jaida, group.vec.Handan))
    Z.block <- list()
    for (i in 1:length(re.block.inds)) Z.block[[i]] <- as.formula(paste("~Z.Jaida[,c(", 
                                                                        paste(re.block.inds[[i]], collapse = ","), ")]-1"))
    if (length(re.block.inds) == 1) {
      if (family == "gaussian")
           lme.fit <- lme(y ~ -1 + X.Declan, random = pdIdent(~-1 + Z.Jaida), data = data.fr, method = spar.method,correlation=correlation,control=control)
      
      if (family != "gaussian") {
        require("MASS")
        offs <- spm.info$off.set
        if(is.null(offs))
        lme.fit <- glmmPQL(y ~ -1 + X.Declan, random = list(group.vec.Handan = pdIdent(~-1 + 
                                                              Z.Jaida)), data = data.fr, family = family, 
                           weights = weights, niter = 30,correlation=correlation,control=control)
        else
            lme.fit <- glmmPQL(y ~ -1 + offset(offs)+X.Declan, random = list(group.vec.Handan = pdIdent(~-1 + 
                                                              Z.Jaida)), data = data.fr, family = family, 
                           weights = weights, niter = 30,correlation=correlation,control=control)
          
      }
    }
    if (length(re.block.inds) > 1) {
      if (family == "gaussian") 
        lme.fit <- lme(y ~ -1 + X.Declan, random = list(group.vec.Handan = pdBlocked(Z.block, 
                                                          pdClass = rep("pdIdent", length(Z.block)))), 
                       data = data.fr, method = spar.method,correlation=correlation,control=control)
      if (family != "gaussian") {
        require("MASS")
        offs <- spm.info$off.set
        if(is.null(offs))
        lme.fit <- glmmPQL(y ~ -1 + X.Declan, random = list(group.vec.Handan = pdBlocked(Z.block, 
                                                              pdClass = rep("pdIdent", length(Z.block)))), 
                           data = data.fr, family = family, weights = weights, 
                           niter = 30,correlation=correlation,control=control)
        else
                  lme.fit <- glmmPQL(y ~ -1 + offset(offs)+X.Declan, random = list(group.vec.Handan = pdBlocked(Z.block, 
                                                              pdClass = rep("pdIdent", length(Z.block)))), 
                           data = data.fr, family = family, weights = weights, 
                           niter = 30,correlation=correlation,control=control)

      }
    }
    lme.fit <- c(lme.fit, list(sigma = summary(lme.fit)$sigma))
  }
  if (is.null(Z)) {
    data.fr <- cbind(y, X.Declan, group.vec.Handan)
    G <- NULL
    if (family == "gaussian") {
      if (!is.null(correlation))
      lm.fit <- gls(y ~ -1 + X.Declan,correlation=correlation)
      else
        lm.fit <- gls(y ~ -1 + X.Declan,correlation=correlation)
      lme.fit <- list(coef = list(fixed = lm.fit$coef), 
                      sigma = summary(lm.fit)$sigma)
    }
    if (family != "gaussian") {
      if (!is.null(spm.info$off.set)) {
        if (!is.null(X)) 
          glm.fit <- glm(y ~ -1 + X.Declan, offset = spm.info$off.set, 
                         family = family)
        if (is.null(X)) 
          glm.fit <- glm(y ~ 1, offset = spm.info$off.set, 
                         family = family)
      }
      if (is.null(spm.info$off.set)) {
        if (!is.null(X)) 
          glm.fit <- glm(y ~ -1 + X.Declan, family = family)
        if (is.null(X)) 
          glm.fit <- glm(y ~ 1, family = family)
      }
      lme.fit <- list()
      lme.fit$coef$fixed <- glm.fit$coef
      lme.fit$coef$random <- NULL
      lme.fit$loglik <- NULL
    }
  }
  RR <- NULL
  if (!is.null(Z)) {
    lme.fit$coef$random <- unlist(lme.fit$coef$random)
    sig.u.hat <- lme.fit$sigma * exp(unlist(lme.fit$modelStruct))
    diag.sqrt.G <- NULL
    for (ib in 1:length(re.block.inds)) diag.sqrt.G <- c(diag.sqrt.G,  rep(sig.u.hat[ib], length(re.block.inds[[ib]])))
    G <- diag(diag.sqrt.G^2)
    }

#update for the correlated errors the covariance matrix

  if(!is.null(Z)&!is.null(correlation))
    RR <- corMatrix(lme.fit$modelStruct$corStruct)
  resid.var <- lme.fit$sigma^2
  if (auto.spar.select == FALSE) {
    if (family == "gaussian") {
      G <- resid.var * diag(diag.G)
      qr.out <- lmeFitQr(y, X, Z, G, resid.var = resid.var)
      coef.ests <- qr.out$coefficients[1:(ncol(X) + ncol(Z))]
      lme.fit <- list()
      lme.fit$coef$fixed <- coef.ests[1:ncol(X)]
      lme.fit$coef$random <- coef.ests[(1 + ncol(X)):length(coef.ests)]
    }
    if ((family != "gaussian") & (!is.null(Z))) {
      G <- diag(diag.G)
      C.mat <- cbind(X, Z)
      ridge.vec <- c(rep(0, ncol(X)), 1/diag.G)
      if (!is.null(spm.info$off.set)) 
        ridge.reg.fit <- airls.ridge(C.mat, y, off.var = spm.info$off.set, 
                                    ridge.vec = ridge.vec, max.it = 50, acc = 1e-06, 
                                    family = family)
      else ridge.reg.fit <- airls.ridge(C.mat, y, ridge.vec = ridge.vec, 
                                       max.it = 50, acc = 1e-06, family = family)
      lme.fit <- list()
      lme.fit$coef$fixed <- ridge.reg.fit$coef[1:ncol(X)]
      lme.fit$coef$random <- ridge.reg.fit$coef[(1 + ncol(X)):ncol(C.mat)]
      lme.fit$loglik <- NULL
    }
  }
  if (auto.spar.select == TRUE) {
    if ((!is.null(spm.info$pen)) | (!is.null(spm.info$krige))) {
      sigu2.hat <- rep(0, length(re.block.inds))
      if (!is.null(Z)) 
        for (ib in 1:length(re.block.inds)) sigu2.hat[ib] <- diag(G)[re.block.inds[[ib]][1]]
      if (is.null(spm.info$krige)) {
        basis.type <- spm.info$pen$basis
        for (ip in 1:ncol(as.matrix(spm.info$pen$x))) {
          deg.val <- spm.info$pen$degree[ip]
          if (basis.type == "trunc.poly") 
            spm.info$pen$spar[[ip]] <- exp(log(resid.var/sigu2.hat[ip])/(2 * 
                                                                         deg.val))
          else spm.info$pen$spar[[ip]] <- exp(log(resid.var/sigu2.hat[ip])/deg.val)
        }
      }
      if (is.null(spm.info$pen)) {
        deg.val <- spm.info$krige$degree
        spm.info$krige$spar <- exp(log(resid.var/sigu2.hat)/deg.val)
      }
      if ((!is.null(spm.info$pen)) & (!is.null(spm.info$krige))) {
        basis.type <- spm.info$pen$basis
        for (ip in 1:ncol(as.matrix(spm.info$pen$x))) {
          deg.val <- spm.info$pen$degree[ip]
          var.rats <- (resid.var/sigu2.hat[ip])
          if (basis.type == "trunc.poly") 
            spm.info$pen$spar[[ip]] <- var.rats^(1/(2 * 
                                                    deg.val))
          else spm.info$pen$spar[[ip]] <- var.rats^(1/(deg.val))
        }
        num.pen <- ncol(as.matrix(spm.info$pen$x))
        var.rat <- (resid.var/sigu2.hat[num.pen + 1])
        deg.val <- spm.info$krige$degree
        spm.info$krige$spar <- exp(log(var.rat)/deg.val)
      }
    }
  }

#if variance of the random effects is given, after
#estimate redefine the basis matrix again to Z
#and update the covariance matrix of the random effects

  if (!is.null(Si.b) & (!is.null(Z))) {
    G <- t(t(G) * Si.b)
    Z <- design.info$Z
    lme.fit$coef$random <- as.vector(unlist(lme.fit$coef$random)) * 
      sqrt(Si.b)
  }

  if (family == "gaussian") 
    aux.info <- almeAux(X, Z, G, RR,resid.var, block.inds)
  if (family != "gaussian") {
    if (is.null(Z)) {
      R <- glm.fit$R
      rinv <- backsolve(R, diag(ncol(X)))
      cov.mat <- rinv %*% t(rinv)
      df.fit <- ncol(X)
      df.res <- length(y) - df.fit
      df <- rep(1, ncol(X))
      random.var <- NULL
      aux.info <- list(cov.mat = cov.mat, df = df, block.inds = block.inds, 
                       random.var = random.var, df.fit = df.fit, df.res = df.res)
    }
    if (!is.null(Z)) {
      if (auto.spar.select == TRUE) {
        C.mat <- cbind(X, Z)
        if(!is.null(RR))
          {
            R.svd <- svd(RR)
            RR12 <- R.svd$u%*%diag(1/sqrt(R.svd$d))%*%t(R.svd$v)
            C.mat <- RR12%*%C.mat
          }
        diag.G <- diag(G)
        ridge.vec <- c(rep(0, ncol(X)), 1/diag.G)
        if (!is.null(spm.info$off.set)) 
          ridge.reg.fit <- airls.ridge(C.mat, y, off.var = spm.info$off.set, 
                                      ridge.vec = ridge.vec, max.it = 50, acc = 1e-06, 
                                      family = family)
        else ridge.reg.fit <- airls.ridge(C.mat, y, ridge.vec = ridge.vec, 
                                         max.it = 50, acc = 1e-06, family = family)
      }
      aux.info <- glmeAux(X, Z, G, block.inds, ridge.reg.fit, 
                          family)
        C.mat <- cbind(X, Z)
    }
  }
  if (!is.null(Z)) {
    coef.ests <- c(lme.fit$coef$fixed, lme.fit$coef$random)
    C.mat <- cbind(X, Z)
  }
  if (is.null(Z)) {
    coef.ests <- lme.fit$coef$fixed
    C.mat <- X
  }
  mins <- NULL
  maxs <- NULL
  for (j in 1:length(block.inds)) {
    fitted.j <- as.matrix(C.mat[, block.inds[[j]]]) %*% as.matrix(coef.ests[block.inds[[j]]])
    mins[j] <- min(fitted.j)
    maxs[j] <- max(fitted.j)
  }
  aux.info <- c(aux.info, list(mins = mins, maxs = maxs))
  if (family == "gaussian") {
    fitted <- as.vector(C.mat %*% coef.ests)
    resids <- y - fitted
    lme.fit$fitted <- fitted
    lme.fit$residuals <- resids
  }
  if (family != "gaussian") {
    eta.hat <- C.mat %*% coef.ests
    mu.hat <- inv.link(eta.hat, family)
    fitted <- mu.hat
    if (family == "binomial") 
      resids <- binomial()$dev.resids(y, mu.hat, rep(1, 
                                                     length(y)))
    if (family == "poisson") 
      resids <- poisson()$dev.resids(y, mu.hat, rep(1, 
                                                    length(y)))
    Dev <- sum(resids^2)
    wt <- dinv.link(eta.hat, family)
    Dev.wls <- sum((y - inv.link(eta.hat, family))^2/wt)
    lme.fit <- c(lme.fit, list(fitted = fitted, residuals = resids, 
                               deviance = Dev, deviance.wls = Dev.wls))
  }
  rm("group.vec.Handan", pos = 1)
  rm("X.Declan", pos = 1)
  if (!is.null(Z)) 
    rm("Z.Jaida", pos = 1)
  names(lme.fit$coef$fixed) <- NULL
  names(lme.fit$coef$random) <- NULL
  spm.fit.obj <- list(fit = lme.fit, info = spm.info, aux = aux.info)
  class(spm.fit.obj) <- "spm"
  return(spm.fit.obj)
}
