eiread <- function (ei.object, ...) 
{
  .betaB <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ei.object$betab
    }
  }
  .betaW <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ei.object$betaw
    }
  }
  .phi <- function (ei.object) 
  {
    ei.object$phi
  }
  .sbetab <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ei.object$sbetab
    }
  }
  
  .sbetaw <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ei.object$sbetaw
    }
  }
    
  .psisims <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ei.object$psi
    }
  }
    
  .bounds <- function (ei.object) 
  {
    bounds1(ei.object$x, ei.object$t, ei.object$n)
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
  
  .CI80b <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      betab <- ei.object$betabs[ok, ]
      lwr <- vector(mode = "numeric", length = length(ei.object$x))
      upr <- vector(mode = "numeric", length = length(ei.object$x))
      lwr[ok] <- apply(betab, 1, function(x) quantile(x, probs = c(0.1)))
      lwr[!ok] <- NA
      upr[ok] <- apply(betab, 1, function(x) quantile(x, probs = c(0.9)))
      upr[!ok] <- NA
      return(cbind(lwr, upr))
    }
  }
    
  .CI80w <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
      betaw <- ei.object$betaws[ok, ]
      lwr <- vector(mode = "numeric", length = length(ei.object$x))
      upr <- vector(mode = "numeric", length = length(ei.object$x))
      lwr[ok] <- apply(betaw, 1, function(x) quantile(x, probs = c(0.1)))
      lwr[!ok] <- NA
      upr[ok] <- apply(betaw, 1, function(x) quantile(x, probs = c(0.9)))
      upr[!ok] <- NA
      return(cbind(lwr, upr))
    }
  }

  .abounds <- function (ei.object) 
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
    
  .aggs <- function (ei.object) 
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
      return(cbind(Bbgg, Bwgg))
    }
  }

  .maggs <- function (ei.object) 
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

  .VCaggs <- function (ei.object) 
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
      vc <- matrix(c(var(Bbgg), cov(Bbgg, Bwgg), cov(Bbgg, 
                                                     Bwgg), var(Bwgg, Bwgg)), nrow = 2)
      return(vc)
    }
  }

  .eaggbias <- function (ei.object) 
  {
    if (!("betabs" %in% names(ei.object))) {
      message("Error: This eiread function requires an ei.sim object.")
    }
    if ("betabs" %in% names(ei.object)) {
      x <- ei.object$x
      mbetab <- ei.object$betab
      mbetaw <- ei.object$betaw
      lm.b <- lm(mbetab ~ x)
      lm.w <- lm(mbetaw ~ x)
      output <- list(summary(lm.b)$coefficients, summary(lm.w)$coefficients)
      names(output) <- c("betaB", "betaW")
      return(output)
    }
  }

  .goodman <- function (ei.object) 
  {
    x <- ei.object$x
    t <- ei.object$t
    lm.g <- lm(t ~ x)
    w <- 1 - x
    lm.w <- lm(t ~ w)
    BetaW <- summary(lm.g)$coefficients[1, ]
    BetaB <- summary(lm.w)$coefficients[1, ]
    coefs <- rbind(BetaB, BetaW)
    return(coefs)
  }
  
  function.list <- list(betab = .betaB, betaw = .betaW, phi = .phi, 
                        sbetab = .sbetab, sbetaw = .sbetaw, psisims = .psisims, 
                        bounds = .bounds, CI80b = .CI80b, CI80w = .CI80w, abounds = .abounds, 
                        aggs = .aggs, maggs = .maggs, VCaggs = .VCaggs, eaggbias = .eaggbias, 
                        goodman = .goodman)
  dec <- ei.object$precision
  arguments <- list(...)
  results <- list()
  for (arg in arguments) {
    if (arg %in% names(function.list)) 
      results[[arg]] <- floor(function.list[[arg]](ei.object) * 
                                10^dec)/10^dec
    else results[[arg]] <- NA
  }
  if (length(results) == 1) 
    results <- results[[1]]
  if (length(results) < 1) 
    warning("qi results object is empty")
  results
}
