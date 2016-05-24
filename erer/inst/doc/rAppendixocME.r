# A new function: marginal effect for ordered choice
ocME <- function(w, rev.dum = TRUE, digits = 3)
{
  # 1. Check inputs; similar to ocProb()
  # 2. Get data out: x may contains factors so use model.matrix
  lev <- w$lev; J <- length(lev)
  x.name <- attr(x = w$terms, which = "term.labels")
  x2 <- w$model[, x.name]
  ww <- paste("~ 1", paste("+", x.name, collapse = " "), collapse = " ")
  x <- model.matrix(as.formula(ww), data = x2)[, -1] # factor x changed
  x.bar <- as.matrix(colMeans(x))
  b.est <- as.matrix(coef(w)); K <- nrow(b.est)
  xb <- t(x.bar) %*% b.est; z <- c(-10^6, w$zeta, 10^6)

  pfun <- switch(w$method, probit = pnorm, logistic = plogis)
  dfun <- switch(w$method, probit = dnorm, logistic = dlogis)
  V2 <- vcov(w)  # increase covarance matrix by 2 fixed thresholds
  V3 <- rbind(cbind(V2, 0, 0), 0, 0)
  ind <- c(1:K, nrow(V3)-1, (K+1):(K+J-1), nrow(V3))
  V4 <- V3[ind,]; V5 <- V4[, ind]

  # 3. Calcualate marginal effects (ME)
  # 3.1 ME value
  f.xb <- dfun(z[1:J] - xb) - dfun(z[2:(J+1)] - xb)
  me <- b.est %*%  matrix(data = f.xb, nrow = 1)
  colnames(me) <- paste("effect", lev, sep = ".")

  # 3.2 ME error
  se <- matrix(0, nrow = K, ncol = J)
  for (j in 1:J) {
    u1 <- c(z[j] - xb); u2 <- c(z[j+1] - xb)
    if (w$method == "probit") {
      s1 <- -u1
      s2 <- -u2
    } else {
      s1 <- 1 - 2 * pfun(u1)
      s2 <- 1 - 2 * pfun(u2)
    }
    d1 <-      dfun(u1) * (diag(1,K,K) - s1 * (b.est %*% t(x.bar)))
    d2 <- -1 * dfun(u2) * (diag(1,K,K) - s2 * (b.est %*% t(x.bar)))
    q1 <-      dfun(u1) * s1 * b.est
    q2 <- -1 * dfun(u2) * s2 * b.est
    dr <- cbind(d1 + d2, q1, q2)
    V <- V5[c(1:K, K+j, K+j+1), c(1:K, K+j, K+j+1)]
    cova <- dr %*% V %*% t(dr)
    se[, j] <- sqrt(diag(cova))
  }
  colnames(se) <- paste("SE", lev, sep = ".")
  rownames(se) <- colnames(x)

  # 4. Revise ME and error for dummy variable
  if (rev.dum) {
    for (k in 1:K) {
      if (identical(sort(unique(x[, k])), c(0, 1))) {
        for (j in 1:J) {
          x.d1 <- x.bar; x.d1[k, 1] <- 1
          x.d0 <- x.bar; x.d0[k, 1] <- 0
          ua1 <-z[j] - t(x.d1) %*% b.est; ub1 <- z[j+1] - t(x.d1) %*% b.est
          ua0 <-z[j] - t(x.d0) %*% b.est; ub0 <- z[j+1] - t(x.d0) %*% b.est
          me[k, j] <- pfun(ub1) - pfun(ua1) - (pfun(ub0) - pfun(ua0))
          d1 <- (dfun(ua1) - dfun(ub1)) %*% t(x.d1) -
                (dfun(ua0) - dfun(ub0)) %*% t(x.d0)
          q1 <- -dfun(ua1) + dfun(ua0); q2 <-  dfun(ub1) - dfun(ub0)
          dr <- cbind(d1, q1, q2)
          V <- V5[c(1:K, K+j, K+j+1), c(1:K, K+j, K+j+1)]
          se[k, j] <- sqrt(c(dr %*% V %*% t(dr)))
        }
      }
    }
  }

  # 5. Output
  t.value <- me / se
  p.value <- 2 * (1 - pt(abs(t.value), w$df.residual))
  out <- list()
  for (j in 1:J) {
    out[[j]] <- round(cbind(effect = me[, j], error = se[, j],
      t.value = t.value[, j], p.value = p.value[, j]), digits)
  }
  out[[J+1]] <- round(me, digits)
  names(out) <- paste("ME", c(lev, "all"), sep = ".")
  result <- listn(w, out)
  class(result) <- "ocME"
  return(result)
}

# Example: The specification is from MASS.
library(erer); library(MASS); data(housing); tail(housing)
reg3 <- polr(formula = Sat ~ Infl + Type + Cont, data = housing, 
  weights = Freq, Hess = TRUE, method = "probit")
m3 <- ocME(w = reg3); m3