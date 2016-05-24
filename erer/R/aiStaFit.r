aiStaFit <- function(y, share, price, expen, shift = NULL, omit = NULL,
  hom = TRUE, sym = TRUE, AR1 = FALSE, rho.sel = c("all", "mean"), ...) {
  # A. Inputs
  if (!inherits(y, "mts")) {stop("class(y) should be 'mts').\n")}
  nShare <- length(share)
  nExoge <- 1 + length(expen) + length(shift)
  nParam <- nExoge + nShare
  nTotal <- (nShare - 1) * nParam

  # B. Restriction for homogeneity and symmetry
  m.h <- matrix(data = 0, nrow = nShare - 1, ncol = nTotal)
  for (i in 1:(nShare - 1)) {
    for (j in 1:nShare) {
      m.h[i, (i - 1) * nParam + nExoge + j] <- 1 
    }
  }  

  m.s <- matrix(0, nrow = (nShare - 2) * (nShare - 1) / 2, ncol = nTotal)
  k <- 0
  for (i in 1:(nShare - 2)) {
    for (j in (i + 1):(nShare - 1)) {
      k <- k + 1
      m.s[k, (i - 1) * nParam + nExoge + j] <- 1
      m.s[k, (j - 1) * nParam + nExoge + i] <- -1
    }
  }

  m.hs <- rbind(m.h, m.s)
  r.s  <- rep(x = 0, times = nrow(m.s))
  r.h  <- rep(x = 0, times = nrow(m.h))
  r.hs <- rep(x = 0, times = nrow(m.hs))
  if (!hom & !sym) {aa <- NULL; bb <- NULL}
  if ( hom & !sym) {aa <- m.h ; bb <- r.h }
  if (!hom &  sym) {aa <- m.s ; bb <- r.s }
  if ( hom &  sym) {aa <- m.hs; bb <- r.hs}

  # C. Estimation
  fa <- bsFormu(name.y = share, name.x = c(shift, expen, price), 
      intercept = TRUE)
  if(is.null(omit)) {
     nOmit <- length(share); omit <- share[nOmit]
  } else { 
    nOmit <- which(share == omit)
  }
  sa <- fa[-nOmit]    
  est <- systemfitAR(formula = sa, data = data.frame(y), method = "SUR",
    restrict.matrix = aa, restrict.rhs = bb, AR1 = AR1, rho.sel = rho.sel)

#   # Omitted equation: coefficient
#    co <- matrix(0, nrow = nShare - 1, ncol = nParam)
#    for (i in 1:(nShare - 1)) {
#      co[i, ] <- coef(est)[((i-1) * nParam + 1):(i * nParam)]
#    }
#    co.omit <- c(1, rep(0, (nParam-1))) - colSums(co)
#   
#   # Omitted equation: variance for gamma_last and exoge variables
#    cof <- coef(est); vco <- vcov(est)
#    gama.vc <- 0
#    for (i in 1:(nShare - 1)) { 
#      for (j in 1:(nShare - 1)) {
#        gama.vc <- gama.vc + vco[nParam * i, nParam * j]
#      }
#    }
#    ex.vc <- rep(0, nExoge)
#    for(k in 1:nExoge) {
#      for (i in 1:(nShare - 1)) { 
#        for (j in 1:(nShare - 1)) {
#          ex.vc[k] <- ex.vc[k] + vco[nParam * (i-1) + k, nParam * (j-1) + k]
#        }
#      }
#    }
#    
#    # Omitted equation: combined results
#    df <- df.residual(est)
#    c.ex <- co.omit[c(1:nExoge, nParam)]
#    e.ex <- sqrt(c(ex.vc, gama.vc))
#    t.ex <- c.ex / e.ex
#    p.ex <- 2 * ( 1- pt(abs(t.ex), df) )
#    ex <- data.frame(cbind(c.ex, e.ex, t.ex, p.ex))
  
  # D. Output
  result <- listn(y, share, price, expen, shift, omit, 
    nOmit, hom, sym, nShare, nExoge, nParam, nTotal, formula = sa, 
    res.matrix = aa, res.rhs = bb, est, AR1, call = sys.call())
  class(result) <- c("aiStaFit", "aiFit")
  return(result)
}

# print method related to aiStaFit, aiDynFit, aiStaHau
print.aiFit <- function(x, ...) {
  if (inherits(x, "aiStaHau")) {
    print(x$ratio)
  } else {
    print(x$est)
  }
}