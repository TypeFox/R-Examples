summary.ei <-function (object, ...) 
{
  
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
  
  if ("psi" %in% names(object)) {
    ei1 <- object
    numb <- ei1$numb
    covs <- as.logical(ifelse(diag(ei1$hessian) == 0, 0, 
                              1))
    sdphi <- sqrt(diag((solve(ei1$hessian[covs, covs]))))
    zbmiss <- ifelse(covs[6:(5 + numb)] == FALSE, TRUE, FALSE)
    zwmiss <- ifelse(covs[(6 + numb):length(covs)] == FALSE, 
                     TRUE, FALSE)
    names <- c("Bb0", "Bw0", "sigB", "sigW", "rho")
    if (zbmiss == FALSE & zwmiss == FALSE) {
      sdphi <- c(sdphi[1:5], 0, 0)
      names <- c(names, "Zb", "Zw")
    }
    if (zbmiss == TRUE & zwmiss == FALSE) {
      sdphi <- c(sdphi[1:5], 0, sdphi[(5 + numb):sum(covs)])
      numw <- length(ei1$phi) - (5 + numb)
      wname <- NULL
      for (i in 1:numw) {
        wname[i] = paste("Zw", (i - 1), sep = "")
      }
      names <- c(names, "Zb0", wname)
    }
    if (zbmiss == FALSE & zwmiss == TRUE) {
      sdphi <- c(sdphi, 0)
      bname <- NULL
      for (i in 1:numb) {
        bname[i] = paste("Zb", (i - 1), sep = "")
      }
      names <- c(names, bname, "Zw0")
    }
    if (zbmiss == TRUE & zwmiss == TRUE) {
      sdphi <- c(sdphi, 0, 0)
      names <- c(names, "Zb0", "Zw0")
    }
    mle <- rbind(ei1$phi, sdphi)
    colnames(mle) <- names
    rownames(mle) <- c("", "")
    n <- length(ei1$betab)
    BB <- mean(ei1$psi[, 1:n])
    BW <- mean(ei1$psi[, (n + 1):(2 * n)])
    SB <- mean(ei1$psi[, ((2 * n) + 1)])
    SW <- mean(ei1$psi[, ((2 * n) + 2)])
    RHO <- mean(ei1$psi[, ((2 * n) + 3)])
    psiu <- t(matrix(c(BB, BW, SB, SW, RHO)))
    colnames(psiu) <- c("BB", "BW", "SB", "SW", "RHO")
    rownames(psiu) <- c("")
    BB <- mean(na.omit(ei1$betabs))
    BW <- mean(na.omit(ei1$betaws))
    SB <- sd(as.vector(na.omit(ei1$betabs)))
    SW <- sd(as.vector(na.omit(ei1$betaws)))
    mat <- na.omit(cbind(as.vector(ei1$betabs), as.vector(ei1$betaws)))
    RHO <- cor(mat[, 1], mat[, 2])
    psit <- t(matrix(c(BB, BW, SB, SW, RHO)))
    colnames(psit) <- c("BB", "BW", "SB", "SW", "RHO")
    rownames(psit) <- c("")
    ab <- matrix(eiabounds(ei1), nrow = 2)
    rownames(ab) <- c("lower", "upper")
    colnames(ab) <- c("betab", "betaw")
    magg <- matrix(eimaggs(ei1), nrow = 2)
    rownames(magg) <- c("Bb", "Bw")
    colnames(magg) <- c("mean", "sd")
    output <- list(ei1$erho, ei1$esigma, ei1$ebeta, n, ei1$resamp, 
                   mle, psiu, psit, ab, magg, ei1$precision)
    names(output) <- c("Erho", "Esigma", "Ebeta", "N", "Resamp", 
                       "Maximum likelihood results in scale of estimation (and se's)", 
                       "Untruncated psi's", "Truncated psi's (ultimate scale)", 
                       "Aggregate Bounds", "Estimates of Aggregate Quantities of Interest", 
                       "precision")
    class(output) <- "summary"
    return(output)
  }
  if (!("psi" %in% names(object))) {
    ei1 <- object
    n <- length(ei1$x)
    numb <- ei1$numb
    covs <- as.logical(ifelse(diag(ei1$hessian) == 0, 0, 
                              1))
    sdphi <- sqrt(diag((solve(ei1$hessian[covs, covs]))))
    zbmiss <- ifelse(covs[6:(5 + numb)] == FALSE, TRUE, FALSE)
    zwmiss <- ifelse(covs[(6 + numb):length(covs)] == FALSE, 
                     TRUE, FALSE)
    names <- c("Bb0", "Bw0", "sigB", "sigW", "rho")
    if (zbmiss == FALSE & zwmiss == FALSE) {
      sdphi <- c(sdphi[1:5], 0, 0)
      names <- c(names, "Zb", "Zw")
    }
    if (zbmiss == TRUE & zwmiss == FALSE) {
      sdphi <- c(sdphi[1:5], 0, sdphi[(5 + numb):sum(covs)])
      numw <- length(ei1$phi) - (5 + numb)
      wname <- NULL
      for (i in 1:numw) {
        wname[i] = paste("Zw", (i - 1), sep = "")
      }
      names <- c(names, "Zb0", wname)
    }
    if (zbmiss == FALSE & zwmiss == TRUE) {
      sdphi <- c(sdphi, 0)
      bname <- NULL
      for (i in 1:numb) {
        bname[i] = paste("Zb", (i - 1), sep = "")
      }
      names <- c(names, bname, "Zw0")
    }
    if (zbmiss == TRUE & zwmiss == TRUE) {
      sdphi <- c(sdphi, 0, 0)
      names <- c(names, "Zb0", "Zw0")
    }
    mle <- rbind(ei1$phi, sdphi)
    colnames(mle) <- names
    rownames(mle) <- c("", "")
    ab <- matrix(eiabounds(ei1), nrow = 2)
    rownames(ab) <- c("lower", "upper")
    colnames(ab) <- c("betab", "betaw")
    output <- list(ei1$erho, ei1$esigma, ei1$ebeta, n, mle, 
                   ab, ei1$precision)
    names(output) <- c("Erho", "Esigma", "Ebeta", "N", "Maximum likelihood results in scale of estimation (and se's)", 
                       "Aggregate Bounds", "precision")
    class(output) <- "summary"
    return(output)
  }
}
