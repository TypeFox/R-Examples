# A new function: probabilities for ordered choice
ocProb <- function(w, nam.c, n = 100, digits = 3)
{
  # 1. Check inputs
  if (!inherits(w, "polr")) {
    stop("Need an ordered choice model from 'polr()'.\n")
  }
  if (w$method != "probit" & w$method != "logistic") {
    stop("Need a probit or logit model.\n")
  }
  if (missing(nam.c)) stop("Need a continous variable name'.\n")

  # 2. Abstract data out
  lev <- w$lev; J <- length(lev)
  x.name <- attr(x = w$terms, which = "term.labels")
  x2 <- w$model[, x.name]
  if (identical(sort(unique(x2[, nam.c])), c(0, 1)) ||
    inherits(x2[, nam.c], what = "factor")) {
    stop("nam.c must be a continuous variable.")
  }

  ww <- paste("~ 1", paste("+", x.name, collapse = " "), collapse = " ")
  x <- model.matrix(as.formula(ww), data = x2)[, -1]
  b.est <- as.matrix(coef(w)); K <- nrow(b.est)
  z <- c(-10^6, w$zeta, 10^6)  # expand it with two extreme thresholds
  z2 <- matrix(data = z, nrow = n, ncol = length(z), byrow = TRUE)

  pfun <- switch(w$method, probit = pnorm, logistic = plogis)
  dfun <- switch(w$method, probit = dnorm, logistic = dlogis)
  V2 <- vcov(w)  # increase covarance matrix by 2 fixed thresholds
  V3 <- rbind(cbind(V2, 0, 0), 0, 0)
  ind <- c(1:K, nrow(V3)-1, (K+1):(K+J-1), nrow(V3))
  V4 <- V3[ind, ]; V5 <- V4[, ind]

  # 3. Construct x matrix and compute xb
  mm <- matrix(data = colMeans(x), ncol = ncol(x), nrow = n, byrow = TRUE)
  colnames(mm) <- colnames(x)
  ran <- range(x[, nam.c])
  mm[, nam.c] <- seq(from = ran[1], to = ran[2], length.out = n)
  xb <- mm %*% b.est
  xb2 <- matrix(data = xb, nrow = n, ncol = J, byrow = FALSE)  # J copy

  # 4. Compute probability by category; vectorized on z2 and xb2
  pp <- pfun(z2[, 2:(J+1)] - xb2) - pfun(z2[, 1:J] - xb2)
  trend <- cbind(mm[, nam.c], pp)
  colnames(trend) <- c(nam.c, paste("p", lev, sep="."))

  # 5. Compute the standard errors
  se <- matrix(data = 0, nrow = n, ncol = J)
  for (i in 1:J) {
    z1 <- z[i] - xb; z2 <- z[i+1] - xb
    d1 <- diag(c(dfun(z1) - dfun(z2)), n, n) %*% mm
    q1 <- - dfun(z1); q2 <-   dfun(z2)
    dr <- cbind(d1, q1, q2)
    V <- V5[c(1:K, K+i, K+i+1), c(1:K, K+i, K+i+1)]
    va <- dr %*% V %*% t(dr)
    se[, i] <- sqrt(diag(va))
  }
  colnames(se) <- paste("Pred_SE", lev, sep = ".")

  # 6. Report results
  t.value <- pp / se
  p.value <- 2 * (1 - pt(abs(t.value), n - K))
  out <- list()
  for (i in 1:J) {
    out[[i]] <- round(cbind(predicted_prob = pp[, i], error = se[, i],
      t.value = t.value[, i], p.value = p.value[, i]), digits)
  }
  out[[J+1]] <- round(x = trend, digits = digits)
  names(out) <- paste("predicted_prob", c(lev, "all"), sep = ".")
  result <- listn(w, nam.c, method=w$method, mean.x=colMeans(x), out, lev)
  class(result) <- "ocProb"; return(result)
}

# Example: include "Freq" to have a continuous variable for demo
library(erer); library(MASS); data(housing); str(housing); tail(housing)
reg2 <- polr(formula = Sat ~ Infl + Type + Cont + Freq, data = housing, 
  Hess = TRUE, method = "probit")
p2 <- ocProb(w = reg2, nam.c = 'Freq', n = 300); p2
plot(p2)