ei_est_gen <- function(cand_vector, race_group, total, rho=10, data,
                       table_names, sample=1000, tomog=F, density_plot=F,...) {
  
  # Package ei functions
  ei <- function (formula, total = NULL, Zb = 1, Zw = 1, id = NA, data = NA, 
                  erho = 0.5, esigma = 0.5, ebeta = 0.5, ealphab = NA, ealphaw = NA, 
                  truth = NA, simulate = TRUE, covariate = NULL, lambda1 = 4, 
                  lambda2 = 2, covariate.prior.list = NULL, tune.list = NULL, 
                  start.list = NULL, sample = 1000, thin = 1, burnin = 1000, 
                  verbose = 0, ret.beta = "r", ret.mcmc = TRUE, usrfun = NULL) 
  {
    dv <- terms.formula(formula)[[2]]
    iv <- terms.formula(formula)[[3]]
    t <- as.character(dv)
    x <- as.character(iv)
    n <- as.character(total)
    id <- as.character(id)
    if (length(dv) == 1) {
      print("Running 2x2 ei")
      if (simulate == FALSE) {
        dbuf <- ei.estimate(t, x, n, id = id, data = data, 
                            Zb = Zb, Zw = Zw, erho = erho, esigma = esigma, 
                            ebeta = ebeta, ealphab = ealphab, ealphaw = ealphaw, 
                            truth = truth)
        return(dbuf)
      }
      if (simulate == TRUE) {
        dbuf <- tryCatch(tryCatch(ei.estimate(t, x, n, id = id, 
                                              data = data, Zb = Zb, Zw = Zw, erho = erho, esigma = esigma, 
                                              ebeta = ebeta, ealphab = ealphab, ealphaw = ealphaw, 
                                              truth = truth), error = function(x) ei(t, x, 
                                                                                     n, id = id, data = data, Zb = Zb, Zw = Zw, erho = 3, 
                                                                                     esigma = esigma, ebeta = ebeta, ealphab = ealphab, 
                                                                                     ealphaw = ealphaw, truth = truth)), error = function(x) ei.estimate(t, 
                                                                                                                                                         x, n, id = id, data = data, Zb = Zb, Zw = Zw, 
                                                                                                                                                         erho = 5, esigma = esigma, ebeta = ebeta, ealphab = ealphab, 
                                                                                                                                                         ealphaw = ealphaw, truth = truth))
        dbuf.sim <- ei.sim(dbuf)
        return(dbuf.sim)
      }
    }
    if (length(dv) > 1) {
      print("Running eiRxC")
      dbuf <- ei.MD.bayes(formula, data = data, total = total, 
                          covariate = covariate, lambda1 = lambda1, lambda2 = lambda2, 
                          covariate.prior.list = covariate.prior.list, tune.list = tune.list, 
                          start.list = start.list, sample = sample, thin = thin, 
                          burnin = burnin, verbose = verbose, ret.beta = ret.beta, 
                          ret.mcmc = ret.mcmc, usrfun = usrfun)
      dbuf$data <- data
      dbuf$total <- n
      dbuf$formula <- formula
      class(dbuf) <- "ei"
      return(dbuf)
    }
  }
  
  # ei.estimate
  ei.estimate <- function (t, x, n, id, Zb = 1, Zw = 1, data = NA, erho = 0.5, 
                           esigma = 0.5, ebeta = 0.5, ealphab = NA, ealphaw = NA, truth = NA, 
                           Rfun = 2, precision = 4) 
  {
    if (!missing(data)) {
      t <- data[[t]]
      x <- data[[x]]
      n <- data[[n]]
      if (is.character(Zb)) 
        Zb <- data[[Zb]]
      if (is.character(Zw)) 
        Zw <- data[[Zw]]
      id <- data[[id]]
    }
    Zb <- as.matrix(Zb)
    Zw <- as.matrix(Zw)
    if (dim(Zb)[1] == 1 & Zb[1, 1] == 1 & dim(Zw)[1] == 1 & Zw[1, 
                                                               1] == 1) 
      Rfun = 5
    if (dim(Zb)[1] == 1 & Zb[1, 1] == 1) 
      Zb <- as.matrix(rep(1, length(x)))
    if (dim(Zw)[1] == 1 & Zw[1, 1] == 1) 
      Zw <- as.matrix(rep(1, length(x)))
    numb <- dim(Zb)[2]
    numw <- dim(Zw)[2]
    start <- c(0, 0, -1.2, -1.2, 0, rep(0, numb + numw))
    message("Maximizing likelihood")
    solution <- ucminf(start, eilike, y = t, x = x, n = n, Zb = Zb, 
                       Zw = Zw, numb = numb, erho = erho, esigma = esigma, ebeta = ebeta, 
                       ealphab = ealphab, ealphaw = ealphaw, Rfun = Rfun, hessian = 3)
    covs <- as.logical(ifelse(diag(solution$hessian) == 0 | diag(solution$hessian) == 
                                1, 0, 1))
    hessian <- solution$hessian[covs, covs]
    output <- list(solution$par, solution$hessian, hessian, erho, 
                   esigma, ebeta, ealphab, ealphaw, numb, x, t, n, Zb, Zw, 
                   truth, precision, covs, Rfun, id)
    names(output) <- c("phi", "hessian", "hessianC", "erho", 
                       "esigma", "ebeta", "ealphab", "ealphaw", "numb", "x", 
                       "t", "n", "Zb", "Zw", "truth", "precision", "covs", "Rfun", 
                       "id")
    class(output) <- "ei"
    return(output)
  }
  
  #ei.sim
  ei.sim <- function (ei.object) 
  {
    hessian <- ei.object$hessianC
    erho <- ei.object$erho
    esigma <- ei.object$esigma
    ebeta <- ei.object$ebeta
    ealphab <- ei.object$ealphab
    ealphaw <- ei.object$ealphaw
    numb <- ei.object$numb
    covs <- ei.object$covs
    Rfun <- ei.object$Rfun
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    Zb <- ei.object$Zb
    Zw <- ei.object$Zw
    truth <- ei.object$truth
    id <- ei.object$id
    precision <- ei.object$precision
    message("Importance Sampling..")
    keep <- matrix(data = NA, ncol = (length(ei.object$phi)))
    resamp <- 0
    while (dim(keep)[1] < 100) {
      keep <- eisamp(t, x, n, Zb, Zw, ei.object$phi, hessian,
                     100, keep, numb = numb, covs, erho, esigma, ebeta, 
                     ealphab, ealphaw, Rfun)
      resamp = resamp + 1
    }
    keep <- keep[2:100, ]
    mu <- keep[, 1:2]
    sd <- keep[, 3:4]
    rho <- keep[, 5]
    Bb0v <- keep[, 6:(5 + numb)]
    Bw0v <- keep[, (6 + numb):length(ei.object$phi)]
    sd[, 1] <- exp(sd[, 1])
    sd[, 2] <- exp(sd[, 2])
    Zb <- as.matrix(Zb)
    Zw <- as.matrix(Zw)
    Bb0v <- as.matrix(Bb0v)
    Bw0v <- as.matrix(Bw0v)
    mu1 <- mu[, 1] * (0.25 + sd[, 1]^2) + 0.5 + t(as.matrix(apply(Zb, 
                                                                  2, function(x) x - mean(x))) %*% t(Bb0v))
    mu2 <- mu[, 2] * (0.25 + sd[, 2]^2) + 0.5 + t(as.matrix(apply(Zw, 
                                                                  2, function(x) x - mean(x))) %*% t(Bw0v))
    rho <- (exp(2 * rho) - 1)/(exp(2 * rho) + 1)
    psi <- cbind(mu1, mu2, sd, rho)
    bb <- psi[, 1:length(x)]
    bw <- psi[, (length(x) + 1):(length(x) * 2)]
    sb <- psi[, (length(x) * 2 + 1)]
    sw <- psi[, (length(x) * 2 + 2)]
    rho <- psi[, (length(x) * 2 + 3)]
    omx <- 1 - x
    sbw <- rho * sb * sw
    betab <- matrix(nrow = length(x), ncol = dim(keep)[1])
    betaw <- matrix(nrow = length(x), ncol = dim(keep)[1])
    homoindx <- ifelse(x == 0, 1, 0)
    homoindx <- ifelse(x == 1, 2, homoindx)
    enumtol = 1e-04
    cT0 <- t < enumtol & homoindx == 0
    cT1 <- t > (1 - enumtol) & homoindx == 0
    ok <- ifelse(homoindx == 0 & cT0 == 0 & cT1 == 0, T, F)
    wh <- homoindx == 1
    bl <- homoindx == 2
    for (i in 1:dim(keep)[1]) {
      sig2 <- sb[i]^2 * x^2 + sw[i]^2 * omx^2 + sbw[i] * 2 * 
        x * omx
      omega <- sb[i]^2 * x + sbw[i] * omx
      eps <- t - (bb[i, ]) * x - (bw[i, ]) * omx
      mbb <- bb[i, ] + omega/sig2 * eps
      vbb <- sb[i]^2 - (omega^2)/sig2
      vbb = ifelse(vbb < 1 * 10^-32, 1e-04, vbb)
      s <- ifelse(vbb >= 0 & vbb != Inf & !is.na(vbb), sqrt(vbb), 
                  NaN)
      bounds <- bounds1(x, t, n)
      out <- NULL
      for (j in 1:length(x[ok])) {
        out[ok][j] <- rtnorm(1, mean = mbb[ok][j], sd = s[ok][j], 
                             lower = bounds[ok, ][j, 1], upper = bounds[ok, 
                                                                        ][j, 2])
      }
      out[wh] <- NA
      out[bl] <- t[bl]
      out[cT1] <- bounds[cT1, 1]
      out[cT0] <- bounds[cT0, 1]
      betab[, i] = out
    }
    omx <- 1 - x
    for (j in 1:length(x[ok])) {
      betabs <- betab[ok, ][j, ]
      betaw[ok, ][j, ] <- t[ok][j]/omx[ok][j] - betabs * x[ok][j]/omx[ok][j]
    }
    if (sum(wh) > 0) {
      betaw[wh, ] <- as.matrix(rep(1, dim(keep)[1])) %*% t(as.matrix(t[wh]))
    }
    if (sum(bl) > 0) {
      betaw[bl, ] <- NA
    }
    if (sum(cT1) > 0) {
      betaw[cT1, ] <- as.matrix(rep(1, dim(keep)[1])) %*% t(as.matrix(bounds[cT1, 
                                                                             3]))
    }
    if (sum(cT0) > 0) {
      betaw[cT0, ] <- as.matrix(rep(1, dim(keep)[1])) %*% t(as.matrix(bounds[cT0, 
                                                                             3]))
    }
    mbetab <- apply(betab, 1, mean)
    mbetaw <- apply(betaw, 1, mean)
    sdbetab <- apply(betab, 1, sd)
    sdbetaw <- apply(betaw, 1, sd)
    output <- list(ei.object$phi, ei.object$hessian, hessian, 
                   psi, mbetab, mbetaw, sdbetab, sdbetaw, betab, betaw, 
                   resamp, erho, esigma, ebeta, ealphab, ealphaw, numb, 
                   x, t, n, Zb, Zw, truth, precision, id)
    names(output) <- c("phi", "hessian", "hessianC", "psi", "betab", 
                       "betaw", "sbetab", "sbetaw", "betabs", "betaws", "resamp", 
                       "erho", "esigma", "ebeta", "ealphab", "ealphaw", "numb", 
                       "x", "t", "n", "Zb", "Zw", "truth", "precision", "id")
    class(output) <- "ei"
    return(output)
  }
  
  #eisamp
  eisamp <- function (t, x, n, Zb, Zw, par, varcv, nsims, keep, numb, covs, 
                      erho, esigma, ebeta, ealphab, ealphaw, Rfun) 
  {
    import1 <- NULL
    varcv2 <- solve(varcv)/4
    draw <- rmvnorm(nsims, par[covs], varcv2)
    varcv3 <- solve(varcv2)
    phiv <- dmvnorm(draw, par[covs], varcv2, log = T)
    zbmiss <- ifelse(covs[6] == FALSE, TRUE, FALSE)
    zwmiss <- ifelse(covs[(6 + numb)] == FALSE, TRUE, FALSE)
    if (zbmiss == TRUE & zwmiss == FALSE) {
      draw <- cbind(draw[, 1:5], rep(1, nsims), draw[, (5 + 
                                                          numb):sum(covs)])
    }
    if (zbmiss == FALSE & zwmiss == TRUE) {
      draw <- cbind(draw, rep(1, nsims))
    }
    if (zbmiss == TRUE & zwmiss == TRUE) {
      draw <- cbind(draw, rep(1, nsims), rep(1, nsims))
    }
    import1 <- apply(as.matrix(1:nsims), 1, function(i) -eilike(as.vector(draw[i, 
                                                                               ]), t, x, n, Zb, Zw, numb = numb, erho, esigma, ebeta, 
                                                                ealphab, ealphaw, Rfun) - phiv[i])
    ok <- !is.nan(import1)
    lnir <- import1 - max(import1[ok])
    ir <- NA
    ir[ok] <- exp(lnir[ok])
    tst <- ifelse(is.finite(ir), ir > runif(1, 0, 1), FALSE)
    keep <- rbind(keep, draw[tst, ])
    return(keep)
  }
  
  #eilike
  eilike <- function (param, y, x, n, Zb, Zw, numb, erho, esigma, ebeta, 
                      ealphab, ealphaw, Rfun) 
  {
    Bb0 <- param[1]
    Bw0 <- param[2]
    sb0 <- param[3]
    sw0 <- param[4]
    rho0 <- param[5]
    Bb0v <- param[6:(5 + numb)]
    Bw0v <- param[(numb + 6):length(param)]
    sb = exp(sb0)
    sw = exp(sw0)
    Zb <- as.matrix(Zb)
    Zw <- as.matrix(Zw)
    bb = Bb0 * (0.25 + sb^2) + 0.5 + as.matrix(apply(Zb, 2, function(x) x - 
                                                       mean(x))) %*% as.matrix(Bb0v)
    bw = Bw0 * (0.25 + sw^2) + 0.5 + as.matrix(apply(Zw, 2, function(x) x - 
                                                       mean(x))) %*% as.matrix(Bw0v)
    rho = (exp(2 * rho0) - 1)/(exp(2 * rho0) + 1)
    sigb2 <- sb^2
    sigw2 <- sw^2
    sigbw = rho * sb * sw
    homoindx <- ifelse(x == 0, 1, 0)
    homoindx <- ifelse(x == 1, 2, homoindx)
    enumtol = 1e-04
    cT0 <- y < enumtol & homoindx == 0
    cT1 <- y > (1 - enumtol) & homoindx == 0
    ok <- ifelse(homoindx == 0 & cT0 == 0 & cT1 == 0, T, F)
    omx <- 1 - x
    mu = bb * x + bw * omx
    epsilon = y - mu
    s2 = sigb2 * (x^2) + sigw2 * (omx^2) + 2 * sigbw * x * omx
    omega = sigb2 * x + sigbw * omx
    ebb = bb + (omega/s2) * epsilon
    vbb = sigb2 - (omega^2)/s2
    vbb = ifelse(vbb < 1 * 10^-32, 1e-04, vbb)
    bounds <- bounds1(x, y, n)
    s <- ifelse(vbb >= 0 & vbb != Inf & !is.na(vbb), sqrt(vbb), 
                NaN)
    res <- NULL
    b.s = (bounds[ok, ][, 2] - ebb[ok])/s[ok]
    as = (bounds[ok, ][, 1] - ebb[ok])/s[ok]
    res[ok] <- log(pnorm(as, lower.tail = F) - pnorm(b.s, lower.tail = F))
    R <- NULL
    bs <- as.matrix(cbind(bb, bw))
    R[ok] <- eicreateR(ok, Rfun, bb, bw, sb, sw, rho, x)
    llik.het <- -0.5 * sum((log(s2[ok]) + (epsilon[ok]^2)/(s2[ok])))
    llik.het <- llik.het + sum(res[ok]) - sum(R[ok])
    wh <- homoindx == 1
    llik.wh = 0
    if (sum(wh) > 0) {
      epsilon = y[wh] - bw[wh]
      llik.whi = -0.5 * (log(sigw2) + (epsilon^2)/(sigw2))
      llik.wh = -0.5 * sum((log(sigw2) + (epsilon^2)/(sigw2)))
      bnds = cbind(rep(0, sum(wh)), rep(1, sum(wh)))
      Ebb = bb[wh] + rho * (sb/sw) * epsilon
      vbb = sigb2 * (1 - rho^2)
      vbb = ifelse(vbb < 1 * 10^-32, 1e-04, vbb)
      s <- ifelse(vbb >= 0 & vbb != Inf & !is.na(vbb), sqrt(vbb), 
                  NaN)
      b.s = (bnds[, 2] - Ebb)/s
      as = (bnds[, 1] - Ebb)/s
      res <- log(pnorm(as, lower.tail = F) - pnorm(b.s, lower.tail = F))
      R[wh] <- eicreateR(wh, Rfun, bb, bw, sb, sw, rho, x)
      llik.wh = llik.wh + sum(res) - sum(R[wh])
    }
    bl <- homoindx == 2
    llik.bl = 0
    if (sum(bl) > 0) {
      epsilon = y[bl] - bb[bl]
      llik.bl = -0.5 * sum((log(sigb2) + (epsilon^2)/(sigb2)))
      bnds = cbind(rep(0, sum(bl)), rep(1, sum(bl)))
      Ebb = bw[bl] + rho * (sw/sb) * epsilon
      vbb = sigw2 * (1 - rho^2)
      vbb = ifelse(vbb < 1 * 10^-32, 1e-04, vbb)
      s <- ifelse(vbb >= 0 & vbb != Inf & !is.na(vbb), sqrt(vbb), 
                  NaN)
      b.s = (bnds[, 2] - Ebb)/s
      as = (bnds[, 1] - Ebb)/s
      res <- log(pnorm(as, lower.tail = F) - pnorm(b.s, lower.tail = F))
      R[bl] <- eicreateR(bl, Rfun, bb, bw, sb, sw, rho, x)
      llik.bl = llik.bl + sum(res) - sum(R[bl])
    }
    llik.cT0 = 0
    if (sum(cT0) > 0) {
      bb.cT0 = bs[cT0, ][1]
      bw.cT0 = bs[cT0, ][2]
      sigma = matrix(c(sigb2, sigbw, sigbw, sigw2), nrow = 2)
      if (sum(cT0) == 1) {
        first = log(dmvnorm(c(0, 0), mean = bs[cT0, ], sigma = sigma))
        second <- eicreateR(cT0, Rfun, bb, bw, sb, sw, rho, 
                            x)
        llik.cT0 = sum(first) - sum(second)
      }
      else {
        first = apply(bs[cT0, ], 1, function(x) log(dmvnorm(c(0, 
                                                              0), mean = as.vector(x), sigma = sigma)))
        second <- NULL
        second <- eicreateR(cT0, Rfun, bb, bw, sb, sw, rho, 
                            x)
        llik.cT0 = sum(first) - sum(second)
      }
    }
    llik.cT1 = 0
    if (sum(cT1) > 0) {
      bb.cT1 = bs[cT1, ][1]
      bw.cT1 = bs[cT1, ][2]
      sigma = matrix(c(sigb2, sigbw, sigbw, sigw2), nrow = 2)
      if (sum(cT1) == 1) {
        first = log(dmvnorm(c(1, 1), mean = bs[cT1, ], sigma = sigma))
        second <- eicreateR(cT1, Rfun, bb, bw, sb, sw, rho, 
                            x)
        llik.cT1 = sum(first) - sum(second)
      }
      if (sum(cT1) > 1) {
        first = apply(as.matrix(bs[cT1, ]), 1, function(x) log(dmvnorm(c(1, 
                                                                         1), mean = as.vector(x), sigma = sigma)))
        second <- NULL
        second <- eicreateR(cT1, Rfun, bb, bw, sb, sw, rho, 
                            x)
        llik.cT1 = sum(first) - sum(second)
      }
    }
    llik = llik.het + llik.bl + llik.wh + llik.cT0 + llik.cT1
    prior = 0
    lpdfnorm = log(dnorm(rho0, 0, sd = erho))
    if (esigma > 0) 
      prior = prior - (1/(2 * esigma^2)) * (sigb2 + sigw2)
    if (erho > 0) 
      prior = prior + lpdfnorm
    if (ebeta > 0 & (mean(bb) < 0)) 
      prior = prior - 0.5 * ((mean(bb)^2)/ebeta)
    if (ebeta > 0 & mean(bb) > 1) 
      prior = prior - 0.5 * ((mean(bb) - 1)^2/ebeta)
    if (ebeta > 0 & mean(bw) < 0) 
      prior = prior - 0.5 * ((mean(bw)^2)/ebeta)
    if (ebeta > 0 & mean(bw) > 1) 
      prior = prior - 0.5 * ((mean(bw) - 1)^2/ebeta)
    if (sum(is.na(ealphab)) == 0) 
      prior = prior + sum(dmvnorm(Bb0v, ealphab[, 1], sigma = diag(ealphab[, 
                                                                           2]^2), log = T))
    if (sum(is.na(ealphaw)) == 0) 
      prior = prior + sum(dmvnorm(Bw0v, ealphaw[, 1], sigma = diag(ealphaw[, 
                                                                           2]), log = T))
    llik = llik + prior
    if (is.na(llik) | abs(llik) == Inf) 
      llik = NaN
    return(-llik)
  }
  
  # eiabounds
  eiabounds <- function (ei.object) 
  {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    bounds <- bounds1(x, t, n)
    omx <- 1 - x
    Nb <- n * x
    Nw <- n * omx
    LAbetaB <- weighted.mean(bounds[, 1], Nb)
    UAbetaB <- weighted.mean(bounds[, 2], Nb)
    LAbetaW <- weighted.mean(bounds[, 3], Nw)
    UAbetaW <- weighted.mean(bounds[, 4], Nw)
    return(matrix(c(LAbetaB, UAbetaB, LAbetaW, UAbetaW), nrow = 2))
  }
  
  # eimaggs
  eimaggs <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      x <- ei.object$x[ok]
      t <- ei.object$t[ok]
      n <- ei.object$n[ok]
      betab <- ei.object$betabs[ok, ]
      betaw <- ei.object$betaws[ok, ]
      omx <- 1 - x
      Nb <- n * x
      Nw <- n * omx
      Bbgg <- vector(mode = "numeric", length = dim(betab)[2])
      for (i in 1:dim(betab)[2]) {
        Bbgg[i] <- weighted.mean(betab[, i], Nb)
      }
      Bwgg <- vector(mode = "numeric", length = dim(betaw)[2])
      for (i in 1:dim(betaw)[2]) {
        Bwgg[i] <- weighted.mean(betaw[, i], Nw)
      }
      return(c(mean(Bbgg), mean(Bwgg), sd(Bbgg), sd(Bwgg)))
    }
  }
  
  # eicreateR
  eicreateR <- function (sub, Rfun, bb, bw, sb, sw, rho, x, numb, numw) 
  {
    out <- NULL
    lower = cbind(-bb[sub]/sb, -bw[sub]/sw)
    upper = cbind(-bb[sub]/sb + 1/sb, -bw[sub]/sw + 1/sw)
    mean = c(0, 0)
    corr = matrix(c(1, rho, rho, 1), nrow = 2)
    if (Rfun == 1) {
      out <- NULL
      makeR <- function(i) {
        qi <- pmvnorm(lower = lower[i, ], upper = upper[i, 
                                                        ], mean = mean, corr = corr)
      }
      out <- foreach(i = 1:length(x[sub]), .combine = "c") %dopar% 
        makeR(i)
      out <- ifelse(out < 0 | out == 0, 1 * 10^-322, out)
      out <- log(out)
      out <- ifelse(is.na(out) | abs(out == Inf), 999, out)
      return(out)
    }
    if (Rfun == 2) {
      makeR <- function(i) {
        qi <- sadmvn(lower = lower[i, ], upper = upper[i, 
                                                       ], mean = mean, varcov = corr)
      }
      out <- apply(as.matrix(1:length(x[sub])), 1, makeR)
      out <- ifelse(out < 0 | out == 0, 1 * 10^-322, out)
      out <- log(out)
      out <- ifelse(is.na(out) | abs(out == Inf), 999, out)
      return(out)
    }
    if (Rfun == 3) {
      fun <- function(x) dmvnorm(x, mean, corr)
      for (i in 1:length(x[sub])) {
        qi <- adaptIntegrate(fun, lowerLimit = lower[i, ], 
                             upperLimit = upper[i, ])$integral
        out[i] <- log(qi)
        out[i] <- ifelse(is.na(out[i]) | abs(out[i] == Inf), 
                         999, out[i])
      }
      return(out)
    }
    if (Rfun == 4) {
      print("Option Rfun==4 is no longer possible")
      stop()
    }
    if (Rfun == 5) {
      lower = lower[1, ]
      upper = upper[1, ]
      qi <- sadmvn(lower = lower, upper = upper, mean = mean, 
                   varcov = corr)
      qi <- ifelse(qi < 1 * 10^-14, 1 * 10^-14, qi)
      qi <- log(qi)
      qi <- ifelse((is.na(qi) | abs(qi) == Inf), 999, qi)
      out <- rep(qi, length(x[sub]))
      return(out)
    }
  }
  # bounds1
  bounds1 <- function (x, t, n) 
  {
    homindx <- NULL
    tx <- NULL
    tomx <- NULL
    LbetaB <- NULL
    UbetaB <- NULL
    LbetaW <- NULL
    UbetaW <- NULL
    omx = 1 - x
    Nb = x * n
    Nw = omx * n
    p = length(x)
    homoindx <- ifelse(x == 0, 1, 0)
    homoindx <- ifelse(x == 1, 2, homoindx)
    tx <- as.matrix(t/x)
    tomx = as.matrix(t/omx)
    tomxx <- as.matrix(tx - (omx/x))
    txx <- as.matrix(tomx - x/(1 - x))
    LbetaB <- apply(tomxx, 1, function(x) max(0, x))
    UbetaB <- apply(tx, 1, function(x) min(x, 1))
    LbetaW <- apply(txx, 1, function(x) max(0, x))
    UbetaW <- apply(tomx, 1, function(x) min(x, 1))
    bl <- homoindx == 2
    LbetaB[bl] = t[bl]
    UbetaB[bl] = t[bl]
    LbetaW[bl] = NA
    UbetaW[bl] = NA
    wh <- homoindx == 1
    LbetaB[wh] = NA
    UbetaB[wh] = NA
    LbetaW[wh] = t[wh]
    UbetaW[wh] = t[wh]
    return(cbind(LbetaB, UbetaB, LbetaW, UbetaW))
  }
  
  # Functions
  list_extract <- function(x) x[,1:2] # sends to lapply to extract indiv column estimates
  
  # Table/Output Row Labeling -- #Added this if statement in to handle calls to one candidate
  seq_split <- 2:length(cand_vector)
  if (length(cand_vector) == 1) {
    rn <- c(cand_vector, "se")
  } else {
    rn <- c(insert(cand_vector, ats= seq_split,values=rep("se",length(cand_vector)-1)), "se")
  }
  # Remove any missing datas
  data <- na.omit(data) 
  
  #Loop Placeholder
  race_group_table <- list()
  
  # Loop over Race Vector
  for (k in 1:length(race_group)) {
    
    # Loop Placeholder
    cand_table <- list() # candidate place holder
    
    # Loop over Candidates
    for (i in 1:length(cand_vector)) {
      
      # Formula object that is looked through
      form <- formula(paste(cand_vector[i], race_group[k])) 
      
      try(ei_out <- ei(form, total = total, erho=rho, data=data, sample=sample,...),silent=T)
      gm <- geterrmessage()
      if(gm == "Maximizing likelihood
         Error in .subset2(x, i, exact = exact) : invalid subscript type 'list'") 
        stop("Maximizing likelihood
             Error in .subset2(x, i, exact = exact) : invalid subscript type 'list'\n
             \n ei package error try re-running ei_est_gen()"
        )
      
      cat(paste("Model:",cand_vector[i], race_group[k], "\n",sep=" "))
      
      print(summary(ei_out))
      
      #Tomography plot
      if (tomog) {
        pdf(paste(cand_vector[i], race_group[k], ".pdf",sep=""))
        plot(ei_out, "tomogE")
        mtext(paste(cand_vector[i], race_group[k], sep=" "), outer=T, line=-1)
        dev.off()
      }
      # Print Out Out density plot for examination
      if(density_plot) {
        pdf(paste("density_plot",k,i,".pdf",sep="_"))
        plot(ei_out, "betab","betaw") # Plot out distribution plots
        mtext(paste(cand_vector[i], race_group[k], sep=" "), outer=T, line=-1)
        dev.off()
      } 
      
      # Extract Beta B and W
      beta_stan_err <- eiread(ei_out, "betab", "sbetab", "betaw", "sbetaw") # beta estimate for minority group
      min_b <- mean(unlist(beta_stan_err[1]), na.rm=T)*100
      min_ste <- mean(unlist(beta_stan_err[2]), na.rm=T)*100
      non_b <- mean(unlist(beta_stan_err[3]), na.rm=T)*100
      non_ste <- mean(unlist(beta_stan_err[4]), na.rm=T)*100
      # Put into useable data frame		
      eimean <- data.frame(c(min_b, min_ste), c(non_b, non_ste)) 
      
      #the results for all candidate are stored here in this list
      cand_table[[i]] <- eimean
      
    } # Close cand_vector loop
    
    cand_table <- rbindlist(cand_table) # cand_table is for one racial group and all candidates
    cand_table <- data.frame(rn, cand_table) # Add in vector for labeling
    
    race_group_table[[k]] <- cand_table # Put candidate results into list
    
  } # Close race group loop
  
  if(length(race_group) == 1) { # For when there is just % Minority vs. % White, for example
    
    race_group_table <- data.frame(race_group_table)
    
  } else{ # For when there are multiple groups (e.g., pct_hisp, pct_asian, pct_white)
    race_group_table <- data.frame( lapply(race_group_table, list_extract) ) # list is length() number of racial groups
    race_group_table <- race_group_table[,c(1,seq(2,ncol(race_group_table),2))] # clean up table
  }
  # Adding on Total Row
  tot <- colSums(race_group_table[seq(1,nrow(race_group_table),2),2:ncol(race_group_table)])
  just_data <- race_group_table[,2:ncol(race_group_table)]
  add <- rbind(just_data, tot)
  add <- data.frame(1:nrow(add), add)
  colnames(add) <- c("Candidate", table_names)
  add[,1] <- c(as.character(race_group_table[,1]), "Total")
  race_group_table <- add 
  
  return(race_group_table)
  
}