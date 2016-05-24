efa <- function(x, ncomp)
{
  nx <- nrow(x)
  Tos <- Fros <- matrix(0, nx, ncomp)
  for (i in 3:nx)
    Tos[i,] <- svd(scale(x[1:i,], scale = FALSE))$d[1:ncomp]
  for (i in (nx-2):1)
    Fros[i,] <- svd(scale(x[i:nx,], scale = FALSE))$d[1:ncomp]

  Combos <- array(c(Tos, Fros[,ncomp:1]), c(nx, ncomp, 2))
  list(forward = Tos, backward = Fros,
       pure.comp = apply(Combos, c(1,2), min))
}

opa <- function(x, ncomp)
{
  Xref <- colMeans(x)
  Xref <- Xref / sqrt(sum(crossprod(Xref))) # scaling

  selected <- rep(0, ncomp)
  for (i in 1:ncomp) {
    Xs <- lapply(1:nrow(x),
                 function(ii, xx, xref) rbind(xref, xx[ii,]),
                 x, Xref)
    dissims <- sapply(Xs, function(xx) det(tcrossprod(xx)))
    selected[i] <- which.max(dissims)
    newX <- x[selected[i],]

    if (i == 1) {
      Xref <- newX / sqrt(crossprod(newX))
    } else {
      Xref <- rbind(Xref, newX / sqrt(sum(crossprod(newX))))
    }
  }
  dimnames(Xref) <- NULL

  list(pure.comp = t(Xref), selected = selected)
}

mcr <- function(x, init, what = c("row", "col"),
                convergence = 1e-8, maxit = 50)
{
  what <- match.arg(what)
  if (what == "col") {
    CX <- init
    SX <- ginv(CX) %*% x
  } else {
    SX <- init
    CX <- x %*% ginv(SX)
  }

  rms <- rep(NA, maxit + 1)
  rms[1] <- sqrt(mean((x - CX %*% SX)^2))

  for (i in 1:maxit) {
    CX <- x %*% ginv(SX)
    SX <- ginv(CX) %*% x

    resids <- x - CX %*% SX
    rms[i+1] <- sqrt(mean(resids^2))
    if ((rms[i] - rms[i+1]) < convergence) break;
  }

  list(C = CX, S = SX, resids = resids, rms = rms[!is.na(rms)])
}
