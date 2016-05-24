roc.area.test <- function(markers, status) {
  if (any(!is.finite(markers))) stop("Marker values should be finite")
  if (any(!is.finite(status))) stop("All status should be finite")
  markers <- as.matrix(markers)[order(status), , drop=FALSE]
  nvar <- ncol(markers)
  nn <- sum(status == 0)
  nd <- sum(status == 1)
  if (min(nn,nd) == 0) stop("Status vector should have least one each of 0 & 1")
  n <- nn + nd
  zzz <- .Fortran("rocarea",
                  as.integer(n),
                  as.integer(nvar),
                  as.integer(nn),
                  as.integer(nd),
                  as.double(markers),
                  area=as.double(rep(0,nvar)),
                  jkarea=as.double(matrix(0,n,nvar)),
                  PACKAGE="clinfun")
  out <- NULL
  out$area <- zzz$area
  areajk <- matrix(zzz$jkarea, n, nvar)
  out$var <- ((nn - 1)^2 * var(areajk[1:nn,  ]))/nn + ((nd - 1)^2 * 
                var(areajk[nn + (1:nd),  ]))/nd
  if (nvar == 2) {
    if (out$area[2] != out$area[1]) {
      out$stat <- (out$area[2] - out$area[1])/sqrt(out$var[1,1] +
         out$var[2,2] - 2*out$var[2,1])
      out$p.value <- 2*pnorm(-abs(out$stat))
    } else {
      out$stat <- 0
      out$p.value <- 1
    }
  }
  if (nvar > 2) {
    A <- diag(1, nvar) - matrix(1, nvar, nvar)/nvar
    x <- (A%*%as.matrix(out$area))[-1]
    v <- (A%*%out$var%*%A)[-1,-1,drop=FALSE]
    if (qr(v)$rank == nvar - 1) {
      out$stat <- sum(solve(v,x)*x)
      out$p.value <- 1 - pchisq(out$stat, df=nvar -1)
      out$df <- nvar - 1
    } else {
      warning("Some markers are perfectly correlated")
      out$stat <- out$p.value <- out$df <- NA
    }
  }
  class(out) <- "roc.area.test"
  out
}

print.roc.area.test <- function(x, ...) {
  if (!inherits(x, 'roc.area.test')) stop("Object not of class roc.area.test")
  k <- length(x$area)
  if (k == 1) {
    cat("  AUC =", x$area, "with s.d", sqrt(x$var), "\n")
  } else {
    if (k == 2) {
      msg <- "from standard normal reference"
    } else {
      msg <- paste("from chi-square (df = ", x$df, ") reference", sep="")
    }
    cat(" ", k, "markers with AUC",  x$area, "\n")
    cat("  test statistic =", x$stat, "\n")
    cat("  p-value =", x$p.value, msg, "\n")
  }
}
