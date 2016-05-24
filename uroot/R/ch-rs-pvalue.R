
# the code for function 'ch.rs.pvalue' below has been 
# ported from gretl code provided by  Ignacio Díaz-Emparanza and M. Paz Moral (2013)
# http://www.ehu.eus/ignacio.diaz-emparanza/packages/Canova_Hansen.gfn
# "Seasonal Stability Tests in gretl. An Application to International Tourism Data"
# Ignacio Díaz-Emparanza and M. Paz Moral, Working Paper Biltoki D.T. 2013.03
# https://addi.ehu.es/handle/10810/10577

ch.rs.pvalue <- function(x, type, lag1, S, n, nobsreg, VMdf)
{
  if (type == "dummy") {
    CMlabel <- if (isTRUE(lag1)) "1_D_qy" else "1_D_q0"
  } else { # type == "trigonometric"  
    tmp <- which(c(1, 2, S-1) == VMdf)
    if (length(tmp) == 0) {
      stop("unexpected value of ", sQuote("VMdf"))
    } else
      CMlabel <- paste0("1_", c("Pi", "Dos", "Ft")[tmp])
      tmp <- if (isTRUE(lag1)) "_qy" else "_q0"
      CMlabel <- paste0(CMlabel, tmp)
  }

  CM <- .CH.CM.tabs[[CMlabel]]

  SdErrs <- CM[,ncol(CM)]
  Coeffs <- CM[,-ncol(CM)]

  rq <- c(0.0001, 0.0002, 0.0005, seq(0.001, 0.01, 0.001), seq(0.015, 0.985, 0.005), 
    seq(0.99, 0.999, 0.001), 0.9995, 0.9998, 0.9999)
  nrq <- length(rq)

  if (type == "dummy") {
    xeplc <- c(1, 1/n, 1/n^2, 1/n^3, S/n, S/n^2, S/n^3)
  } else { # type == "trigonometric"  
    if (VMdf %in% c(1, 2)) {
      xeplc <- c(1, 1/n, 1/n^2, S/n, S/n^2)
    } else # VMdf == S-1, all seasonal frequencies
      xeplc <- c(1, 1/n, 1/n^2, S/n, S/n^2, S-1, (S-1)^2)
  }

  Q1 <- sort(Coeffs %*% xeplc)

  if (x < Q1[1]) {
    return(1)
  } else if (x > Q1[nrq])
    return(0)

  Mdif <- cbind(seq_len(nrq), abs(Q1 - x))
  Msort <- Mdif[order(Mdif[,2]),]
  mascer <- Msort[1,1]
  swindow <- (nobsreg - 1)/2
  minq <- if (mascer <= swindow) 1 else mascer - swindow
  maxq <- if (mascer + swindow >= nrq) nrq else mascer + swindow
  maxq <- if (maxq == 1) nobsreg else maxq
  minq <- if (maxq == nrq) maxq - nobsreg + 1 else minq

  qi <- Q1[seq.int(minq, maxq)]
  pri <- rq[seq.int(minq, maxq)]
  si <- SdErrs[seq.int(minq, maxq)]

  Y <- qchisq(pri, df = 2)
  X <- cbind(1, qi, qi^2, qi^3)

  pri1 <- sqrt(pri / (1 - pri))
  pri2 <- si / pri1
#  Sigma <- matrix(0, nrow = nobsreg, ncol = nobsreg)
  Sigma <- matrix(0, nrow = length(si), ncol = length(si))
  uid <- upper.tri(Sigma, diag = TRUE)
  Sigma[uid] <- tcrossprod(pri1 * si, pri2)[uid]
  Sigma <- t(Sigma) + Sigma - diag(si^2)

  Pinv <- t(solve(chol(Sigma)))
  PY <- Pinv %*% Y
  PX <- Pinv %*% X

  fit <- lm(PY ~ 0 + PX)
  valorcomp <- abs(sum(coef(fit) * c(1, x, x^2, x^3)))

  pchisq(q = valorcomp, df = 2, lower.tail = FALSE)
}
