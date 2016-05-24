ddelap <-
  function (x, alpha, beta, lambda, log = FALSE) {
    DDLAP <- double(length(x))
    DDLAP <- ddelap_C(x, alpha, beta, lambda, log)
    return(DDLAP)
  }

pdelap <-
  function (q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE) {
    PDLAP <- double(length(q))
    PDLAP <- pdelap_C(q, alpha, beta, lambda, lower.tail, log.p)
    return(PDLAP)
  }

qdelap <-
  function (p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE, exact = TRUE) {
    QDLAP <- double(length(p))
    if (exact) {
      QDLAP <- qdelap_C(p, alpha, beta, lambda, lower.tail, log.p)
    } else {
      if (log.p) p <- exp(p)
      if (!lower.tail) p <- 1 - p
      pValid <- p[p > 0 & p < 1]
      pNan <- p[p < 0]
      p0 <- p[p == 0]
      pInf <- p[p >= 1]
      n <- min(10 ^ (ceiling(log(alpha * beta + lambda, 10)) + 5), 1e7)
      NB <- rnbinom(n, mu = alpha * beta, size = alpha)
      P <- rpois(n, lambda = lambda)
      DP <- NB + P
      QValid <- as.vector(quantile(DP, pValid, na.rm = TRUE, type = 8))
      QNan <- rep.int(NaN, times = length(pNan))
      Q0 <- rep.int(0, times = length(p0))
      QInf <- rep.int(Inf, times = length(pInf))
      QDLAP <- as.vector(c(QNan, Q0, QValid, QInf), mode = 'integer')
    }  
    return(QDLAP)
  }

rdelap <-
  function (n, alpha, beta, lambda, exact = TRUE) {
    RDLAP <- double(length(n))
    if (exact) {
      RDLAP <- rdelap_C(n, alpha, beta, lambda)
    } else {
      NB <- rnbinom(max(1e7, n), mu = alpha * beta, size = alpha)
      P <- rpois(max(1e7, n), lambda = lambda)
      DP <- NB + P
      if (n >= 1e7) {
        RDLAP <- DP
      } else {
        RDLAP <- sample(x = DP, size = n, replace = TRUE)
      }
    }
    return(RDLAP)
  }

MoMdelap <- function (x) {
    MoMDLAP <- double(3)
    MoMDLAP <- MoMdelap_C(x)
    if (any(MoMDLAP < 0)) stop ("Method of moments not appropriate for this data; results include negative parameters.")
    return(MoMDLAP)
  }