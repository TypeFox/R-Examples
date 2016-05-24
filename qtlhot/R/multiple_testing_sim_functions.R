MyMultipleTestingSim <- function(nSim,
                          n.ind,
                          n.pheno,
                          add.eff.range = c(0, 0),
                          dom.eff.range = c(0, 0),
                          beta.range = c(0, 0),
                          sig2.1.range = c(1, 1),
                          sig2.2.range = c(1, 1),
                          eq.spacing = FALSE,
                          cross.type = "f2",
                          thr,
                          normalize = FALSE) {
  n.targets <- rep(NA, nSim)
  R2s <- vector(mode = "list", length = nSim)
  BICs <- vector(mode = "list", length = nSim)
  AICs <- vector(mode = "list", length = nSim)
  pvals.j.BIC <- vector(mode = "list", length = nSim)
  pvals.p.BIC <- vector(mode = "list", length = nSim)
  pvals.np.BIC <- vector(mode = "list", length = nSim)
  pvals.j.AIC <- vector(mode = "list", length = nSim)
  pvals.p.AIC <- vector(mode = "list", length = nSim)
  pvals.np.AIC <- vector(mode = "list", length = nSim)
  k <- 1
  while (k <= nSim) { 
    beta <- runif(1, beta.range[1], beta.range[2])
    add.eff <- runif(1, add.eff.range[1], add.eff.range[2])
    dom.eff <- runif(1, dom.eff.range[1], dom.eff.range[2])
    sig2.1 <- runif(1, sig2.1.range[1], sig2.1.range[2])
    sig2.2 <- runif(1, sig2.2.range[1], sig2.2.range[2])
    Cross <- SimCrossCausal(n.ind, rep(100, 3),  101,
                            rep(beta, n.pheno + 1), add.eff, dom.eff, 
                            sig2.1, sig2.2, eq.spacing, cross.type, 
                            normalize)
    Cross <- calc.genoprob(Cross, step = 2)
    scan <- scanone(Cross, pheno.col = 1:(n.pheno + 1), method = "hk")
    su <- summary(scan[, 1:3], thr = thr)
    if (length(su[[1]]) > 0) {
      cat("sim ", k, "\n")
      targets <- which(scan[row.names(su)[1], -c(1:2)] >= thr)[-1]
      n.targets[k] <- length(targets)
      cat("n targets =", n.targets[k], "\n")
      aux <- try(CMSTtestsList(cross = Cross, 
                           pheno1 = "y1", 
                           pheno2 = paste("y", targets, sep = ""),
                           Q.chr = su[[1]][1],
                           Q.pos = su[[2]][1],
                           method = "all",
                           penalty = "both"), silent = TRUE)
      if (class(aux) != "try-error") {
        R2s[[k]] <- aux$R2s
        BICs[[k]] <- aux$BIC.stats
        AICs[[k]] <- aux$AIC.stats
        pvals.j.BIC[[k]] <- aux$pvals.j.BIC
        pvals.p.BIC[[k]] <- aux$pvals.p.BIC
        pvals.np.BIC[[k]] <- aux$pvals.np.BIC
        pvals.j.AIC[[k]] <- aux$pvals.j.AIC
        pvals.p.AIC[[k]] <- aux$pvals.p.AIC
        pvals.np.AIC[[k]] <- aux$pvals.np.AIC
        k <- k + 1
      }
    }
  }
  list(n.targets = n.targets,
       R2s.list = R2s,
       BICs.list = BICs,
       AICs.list = AICs,
       pvals.j.BIC.list = pvals.j.BIC,
       pvals.p.BIC.list = pvals.p.BIC,
       pvals.np.BIC.list = pvals.np.BIC,
       pvals.j.AIC.list = pvals.j.AIC,
       pvals.p.AIC.list = pvals.p.AIC,
       pvals.np.AIC.list = pvals.np.AIC)
}





MyMultipleTestingSim2 <- function(nSim,
                          n.ind,
                          n.pheno,
                          add.eff.1.range = c(0, 0),
                          dom.eff.1.range = c(0, 0),
                          add.eff.h.range = c(0, 0),
                          dom.eff.h.range = c(0, 0),
                          beta.range = c(0, 0),
                          sig2.1.range = c(1, 1),
                          sig2.2.range = c(1, 1),
                          sig2.h.range = c(1, 1),
                          eq.spacing = TRUE,
                          cross.type = "f2",
                          thr,
                          normalize = FALSE) {
  n.targets <- rep(NA, nSim)
  R2s <- vector(mode = "list", length = nSim)
  BICs <- vector(mode = "list", length = nSim)
  AICs <- vector(mode = "list", length = nSim)
  pvals.j.BIC <- vector(mode = "list", length = nSim)
  pvals.p.BIC <- vector(mode = "list", length = nSim)
  pvals.np.BIC <- vector(mode = "list", length = nSim)
  pvals.j.AIC <- vector(mode = "list", length = nSim)
  pvals.p.AIC <- vector(mode = "list", length = nSim)
  pvals.np.AIC <- vector(mode = "list", length = nSim)
  k <- 1
  while (k <= nSim) { 
    beta <- runif(1, beta.range[1], beta.range[2])
    add.eff.1 <- runif(1, add.eff.1.range[1], add.eff.1.range[2])
    dom.eff.1 <- runif(1, dom.eff.1.range[1], dom.eff.1.range[2])
    add.eff.h <- runif(1, add.eff.h.range[1], add.eff.h.range[2])
    dom.eff.h <- runif(1, dom.eff.h.range[1], dom.eff.h.range[2])
    sig2.1 <- runif(1, sig2.1.range[1], sig2.1.range[2])
    sig2.2 <- runif(1, sig2.2.range[1], sig2.2.range[2])
    sig2.h <- runif(1, sig2.h.range[1], sig2.h.range[2])
    Cross <- SimCrossIndep(n.ind, rep(100, 3),  101, rep(beta, n.pheno + 1), 
                           add.eff.1, dom.eff.1, add.eff.h, dom.eff.h, 
                           sig2.1, sig2.2, sig2.h, eq.spacing, cross.type, 
                           normalize)
    Cross <- calc.genoprob(Cross, step = 0)
    scan <- scanone(Cross, pheno.col = 1:(n.pheno + 1), method = "hk")
    su <- summary(scan[, 1:3], thr = thr)
    if (length(su[[1]]) > 0) {
      cat("sim ", k, "\n")
      targets <- which(scan[row.names(su)[1], -c(1:2)] >= thr)[-1]
      n.targets[k] <- length(targets)
      cat("n targets =", n.targets[k], "\n")
      aux <- try(CMSTtestsList(cross = Cross, 
                           pheno1 = "y1", 
                           pheno2 = paste("y", targets, sep = ""),
                           Q.chr = su[[1]][1],
                           Q.pos = su[[2]][1],
                           method = "all",
                           penalty = "both"), silent = TRUE)
      if (class(aux) != "try-error") {
        R2s[[k]] <- aux$R2s
        BICs[[k]] <- aux$BIC.stats
        AICs[[k]] <- aux$AIC.stats
        pvals.j.BIC[[k]] <- aux$pvals.j.BIC
        pvals.p.BIC[[k]] <- aux$pvals.p.BIC
        pvals.np.BIC[[k]] <- aux$pvals.np.BIC
        pvals.j.AIC[[k]] <- aux$pvals.j.AIC
        pvals.p.AIC[[k]] <- aux$pvals.p.AIC
        pvals.np.AIC[[k]] <- aux$pvals.np.AIC
        k <- k + 1
      }
    }
  }
  list(n.targets = n.targets,
       R2s.list = R2s,
       BICs.list = BICs,
       AICs.list = AICs,
       pvals.j.BIC.list = pvals.j.BIC,
       pvals.p.BIC.list = pvals.p.BIC,
       pvals.np.BIC.list = pvals.np.BIC,
       pvals.j.AIC.list = pvals.j.AIC,
       pvals.p.AIC.list = pvals.p.AIC,
       pvals.np.AIC.list = pvals.np.AIC)
}

SimCrossIndep <- function(n.ind, len, n.mar, beta, add.eff.1, dom.eff.1,
                          add.eff.h, dom.eff.h, sig2.1 = 1, sig2.2 = 1, sig2.h = 1, 
                          eq.spacing = FALSE, cross.type = "f2", 
                          normalize = FALSE) {
  n.traits <- length(beta)
  beta <- matrix(rep(beta, each = n.ind), n.ind, n.traits)
  Map <- sim.map(len, n.mar, eq.spacing = eq.spacing, include.x = FALSE)
  Cross <- sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- pull.geno(Cross)
  q1 <- mygeno[, "D1M51"]
  q2 <- mygeno[, "D1M52"]
  if (cross.type == "bc") {
    add.q1 <- q1 - 1.5
    add.q2 <- q2 - 1.5
    y1 <- add.q1 * add.eff.1 + rnorm(n.ind, 0, sqrt(sig2.1))
    h <- add.q2 * add.eff.h + rnorm(n.ind, 0, sqrt(sig2.h))
  }
  if (cross.type == "f2") {
    add.q1 <- q1 - 2
    dom.q1 <- (1 + add.q1) * (1 - add.q1) - 0.5
    y1 <- add.q1 * add.eff.1 + dom.q1 * dom.eff.1 + rnorm(n.ind, 0, sqrt(sig2.1))
    add.q2 <- q2 - 2
    dom.q2 <- (1 + add.q2) * (1 - add.q2) - 0.5
    h <- add.q2 * add.eff.h + dom.q2 * dom.eff.h + rnorm(n.ind, 0, sqrt(sig2.h))
  }
  y <- beta * h + matrix(rnorm(n.ind * n.traits, 0, sqrt(sig2.2)), n.ind, n.traits)
  y <- data.frame(y1, y)
  names(y) <- paste("y", 1 : (n.traits + 1), sep = "")
  if (normalize) {
    apply(y, 2, normal.trans)
  }
  Cross$pheno <- y
  Cross
}





GetPowerFDR <- function(dat, alpha, true.model = "M1") {
  GetCounts <- function(M, alpha) {
    M1 <- sum(M[,1] <= alpha & M[,2] > alpha & M[,3] > alpha & M[,4] > alpha)
    M2 <- sum(M[,1] > alpha & M[,2] <= alpha & M[,3] > alpha & M[,4] > alpha)
    M3 <- sum(M[,1] > alpha & M[,2] > alpha & M[,3] <= alpha & M[,4] > alpha)
    M4 <- sum(M[,1] > alpha & M[,2] > alpha & M[,3] > alpha & M[,4] <= alpha)
    no.call <- nrow(M) - M1 - M2 - M3 - M4
    output <- c(M1, M2, M3, M4, no.call)
  }
  le <- length(dat)
  counts <- matrix(NA, le, 5)
  for (i in 1:le) {
    counts[i, ] <- GetCounts(dat[[i]], alpha)
  }
  counts.sum <- apply(counts, 2, sum)
  calls <- sum(counts.sum[1:4])
  all <- sum(counts.sum[1:5])
  if (true.model == "M1") {
    Power <- sum(counts.sum[1])/all
    FDR <- 0
    if (calls > 0) {
      FDR <- sum(counts.sum[2:4])/calls
    }
  }
  if (true.model == "M3") {
    Power <- sum(counts.sum[3])/all
    FDR <- 0
    if (calls > 0) {
      FDR <- sum(counts.sum[c(1, 2, 4)])/calls
    }
  }
  list(FDR = FDR, Power = Power, counts.sum = counts.sum, counts = counts)
}




GetPowerFdrMatrices <- function(dat.list, alphas, mnms = c("j.bic", "p.bic", 
                         "np.bic", "j.aic", "p.aic", "np.aic"), 
                         true.model = "M1") {
  ll <- length(dat.list)
  la <- length(alphas)
  Power <- FDR <- matrix(NA, ll, la)
  dimnames(FDR) <- list(mnms, as.character(alphas))
  dimnames(Power) <- list(mnms, as.character(alphas))
  for (i in 1:ll) {
  cat("method", i, "\n")
    for (j in 1:la) {
      aux <- GetPowerFDR(dat.list[[i]], alphas[j], true.model)
      FDR[i, j] <- aux[[1]]
      Power[i, j] <- aux[[2]]
      cat("alpha", j, "\n")
    }
  }
  list(FDR = FDR, Power = Power)
}




GetPvalsDistr <- function(dat, n.targets) {
  n <- length(n.targets)
  pvals.1 <- pvals.2 <- pvals.3 <- pvals.4 <- 
    rep(NA, sum(n.targets))
  i.start <- 1
  for (k in 1:n) {
      i.end <- n.targets[k] + i.start - 1
      pvals.1[i.start:i.end] <- dat[[k]][, 1]
      pvals.2[i.start:i.end] <- dat[[k]][, 2]
      pvals.3[i.start:i.end] <- dat[[k]][, 3]
      pvals.4[i.start:i.end] <- dat[[k]][, 4]
      i.start <- i.end + 1
  }
  list(pvals.1 = pvals.1, 
       pvals.2 = pvals.2,
       pvals.3 = pvals.3,
       pvals.4 = pvals.4) 
}






FDRplot <- function(r2s.list, 
                    pvals.list = NULL,
                    scores.list = NULL,
                    xlim = c(0, 1),
                    ylim = c(0, 1), 
                    alpha = 0.05, 
                    main = "", 
                    cex = 0.1,
                    cex.axis = 1.5,
                    cex.lab = 1.5,
                    cex.main = 1.5) {
  xaxis <- seq(xlim[1], xlim[2], length.out = 100)
  yaxis <- seq(ylim[1], ylim[2], length.out = 100)
  par(mar = c(5, 6, 4, 2) + 0.1)
  plot(xaxis, yaxis, type = "n", cex = 0.1, xlim = xlim, ylim = ylim, 
     main = main, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
     xlab = expression(paste(R^{2}*(list(Y[1],Q)))),
     ylab = expression(paste(R^{2}*(list(Y[k],Q)))))
  abline(a=0, b = 1, col = "grey", lwd = 3)
  if (!is.null(pvals.list) & is.null(scores.list)) {
    le1 <- length(pvals.list)
    for (i in 1:le1) {
      cat("", i, "\n")
      le2 <- nrow(pvals.list[[i]])
      for (j in 1:le2) {
        if (pvals.list[[i]][j, 1] <= alpha) {
          points(r2s.list[[i]][j, 1], r2s.list[[i]][j, 2], cex = cex, 
                 col = "blue")
        }
        if (pvals.list[[i]][j, 2] <= alpha) {
          points(r2s.list[[i]][j, 1], r2s.list[[i]][j, 2], cex = cex, 
                 col = "red")
        }
        if (pvals.list[[i]][j, 3] <= alpha) {
          points(r2s.list[[i]][j, 1], r2s.list[[i]][j, 2], cex = cex, 
                 col = "green")
        }
        if (pvals.list[[i]][j, 4] <= alpha) {
          points(r2s.list[[i]][j, 1], r2s.list[[i]][j, 2], cex = cex, 
                 col = "black")
        }
      }
    }
  }
  if (is.null(pvals.list) & !is.null(scores.list)) {
    le1 <- length(scores.list)
    for (i in 1:le1) {
      cat("", i, "\n")
      le2 <- nrow(scores.list[[i]])
      for (j in 1:le2) {
        aux <- which.min(scores.list[[i]][j, 1:4])
        if (aux == 1) {
          points(r2s.list[[i]][j, 1], r2s.list[[i]][j, 2], cex = cex, 
                 col = "blue")
        }
        if (aux == 2) {
          points(r2s.list[[i]][j, 1], r2s.list[[i]][j, 2], cex = cex, 
                 col = "red")
        }
        if (aux == 3) {
          points(r2s.list[[i]][j, 1], r2s.list[[i]][j, 2], cex = cex, 
                 col = "green")
        }
        if (aux == 4) {
          points(r2s.list[[i]][j, 1], r2s.list[[i]][j, 2], cex = cex, 
                 col = "black")
        }
      }
    }
  }
}





GetSubsetForPlot <- function(n, N = 1000, r2.list) {
  aux1 <- rep(NA, N)
  for (i in 1:N) {
    aux1[i] <- r2.list[[i]][1, 1]
  }
  aux2 <- order(aux1)
  NN <- n * round(N/n)
  aux3 <- seq(1, NN, by = round(N/n))
  aux3[n] <- N
  index <- rep(NA, n)
  for (i in 1:n) {
    index[i] <- which(aux2 == aux3[i])
  }
  list(index = index, r2s = aux1[index])
}







