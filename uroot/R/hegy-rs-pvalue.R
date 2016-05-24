
# the code for function 'hegy.rs.pvalue' below has been 
# ported from gretl code provided by Ignacio Díaz-Emparanza (2014)
# http://www.ehu.eus/ignacio.diaz-emparanza/packages/GHegy.gfn
# “Numerical Distribution Functions for Seasonal Unit Root Tests”, 
# Computational Statistics and Data Analysis, Volume 76, pages 237-247 
# The Annals of Computational and Financial Econometrics - 2nd Issue
# DOI: http://dx.doi.org/10.1016/j.csda.2013.03.006

hegy.rs.pvalue <- function(x, type = c("zero", "pi", "pair", "seasall", "all"), 
  deterministic = c(1,0,0), lag.method = c("fixed", "AIC", "BIC"), #"HQC"
  lag.order, S, n, nobsreg)
{
  type <- match.arg(type)
  if (is.numeric(deterministic))
    deterministic <- paste(deterministic, collapse = "")
  lag.method <- match.arg(lag.method)

  if (deterministic == "000")
    stop("the combination of deterministic components ", sQuote(deterministic), 
      " is not available for ", sQuote("pvalue=\"RS\""))

  ##NOTE 
  #casedet = "000" is not considered

  CMlabel <- switch(type, "zero" = "Ct1", "pi" = "Ct2", "pair" = "CF", "seasall" = "CFs", "all" = "CFt")  
  tmp <- switch(deterministic, "100" = "c", "110" = "ct", "101" = "cD", "111" = "cDt")
  CMlabel <- paste(CMlabel, tmp, sep = "_")
  tmp <- switch(lag.method, "fixed" = "fijo", "AIC" = "AIC", "BIC" = "BIC", "HQC" = "HQC")
  CMlabel <- paste(CMlabel, tmp, sep = "_")

  C <- .HEGY.CM.tabs[[CMlabel]]

  nc <- ncol(C)
  C1 <- C[,seq_len(nc-1)]
  sdC1 <- C[,nc]

  r1 <- c(0.0001, 0.0002, 0.0005, seq_len(10)/1000)
  rq <- c(r1, seq.int(15, 985, 5)/1000, 1 - rev(r1))

  xeplc <- c(1, 1/n, 1/n^2, 1/n^3, lag.order/n, lag.order/n^2, lag.order/n^3,
    lag.order^2/n, lag.order^2/n^2, lag.order^2/n^3,
    lag.order^3/n, lag.order^3/n^2, lag.order^3/n^3, S/n, S/n^2, S/n^3)

  Q1 <- C1 %*% xeplc

  nrq <- length(rq)
  centro <- floor(nobsreg/2) + 1

  isFtest <- !(type %in% c("zero", "pi"))

  if (isFtest) {
    stopifnot(NCOL(Q1) == 1)
    Q1 <- rev(Q1)
  }

  bin <- rep(0, nrq)
  bin[Q1 >= 0] <- 1
  sbin <- sum(bin)

  if (isFtest)
    sdC1 <- rev(sdC1)

  Q1s <- sort(Q1)

  if (x < Q1s[1]) {
    if (isFtest) return(1) else return(0)
  } else if (x > tail(Q1s, 1)) {
    if (isFtest) return(0) else return(1)
  } else 
  {
    masque <- which(x > Q1)
    masque <- if (length(masque) > 0) max(masque) else 1
    if (masque < nrq) {
      mascer <- if ((x - Q1[masque]) < (Q1[masque+1] - x)) masque else masque + 1
    } else mascer <- nrq

    centroup <- nrq - centro + 1
    mascer <- if (mascer <= centro) centro else mascer
    mascer <- if (mascer >= centroup) centroup else mascer

    # local regressions with "nobsreg" observations
    qi <- pri <- si <- rep(NA, nobsreg)
    for (i in seq_len(nobsreg))
    {
      j <- mascer - centro + i
      qi[i] <- Q1[j]
      pri[i] <- rq[j]
      si[i] <- sdC1[j]
    }

    Y <- if (isFtest) qchisq(pri, df = 2) else qnorm(pri)
    X <- cbind(1, qi, qi^2, qi^3)
    co <- lm.fit(X, Y)$coef

    # substitute by FGLS
    Sigma <- matrix(0, nrow = nobsreg, ncol = nobsreg)
    for (i in seq_len(nobsreg))
    {
      for (j in seq_len(nobsreg))
        if (i <= j) {
          Sigma[i,j] <- (si[i]*si[j])*sqrt((pri[i]*(1-pri[j]))/(pri[j]*(1-pri[i])))
        } else 
          Sigma[i,j] <- (si[i]*si[j])*sqrt((pri[j]*(1-pri[i]))/(pri[i]*(1-pri[j])))
    }

    Pinv <- t(solve(chol(Sigma)))
    PY <- Pinv %*% Y
    PX <- Pinv %*% X
    co <- lm.fit(PX, PY)$coef

    valorcomp <- sum(co * c(1, x, x^2, x^3))
    if (isFtest) {
      pval <- pchisq(q = abs(valorcomp), df = 2, lower.tail = FALSE)
    } else
      pval <- pnorm(q = valorcomp, lower.tail = TRUE)
  }

  pval
}
