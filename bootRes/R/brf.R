brf <- function(g, p, sb, vnames, ci = 0.05) {
  n <- length(g)
  m <- dim(p)[2]
  param.matrix <- matrix(NA, nrow = m, ncol = 1000)
  if (sb) { # initialize status bar (if TRUE)
    pb <- txtProgressBar(min = 1,  max = 1000, style = 3)
  } 
  for (i in 1:1000) {
    boot.sample <- sample(1:n, n, replace = TRUE)
    boot.g <- g[boot.sample]
    boot.p <- p[boot.sample, ]
    boot.g <- (boot.g - mean(boot.g))/sd(boot.g) # standardize
    boot.p <- apply(boot.p, 2, function(x) { (x - mean(x))/sd(x) }) # standardize
    cor.mat <- cor(boot.p) # correlation matrix X'X (q*q)
    eigen.decomp <- eigen(cor.mat) # eigenvector decomposition
    eigenvectors <- eigen.decomp$vectors # normalized eigenvectors
    eigenvalues <- eigen.decomp$values
    cumprods <- cumprod(eigenvalues) # PVP criterion: calculate cumulative eigenvalues until value < 1
    reduced.eigenvectors <- eigenvectors[, cumprods > 1] # matrix of reduced eigenvectors (q*m)
    pc.scores <- boot.p %*% reduced.eigenvectors # calculate princ comp scores (n*m)
    k <- qr.solve(pc.scores, boot.g) # calculate solution for Z*K = Y (coefficients) (m*1)
    zeros <- rep(0, length(which(cumprods < 1))) # pad K with zero so that Kq*1
    k <- c(k, zeros) # (q*1)
    b <- eigenvectors %*% k # response coefficients (q*1)
    param.matrix[, i] <- b
    if (sb) # update status bar (if TRUE)
      setTxtProgressBar(pb, i)
  }
  brf.coef <- apply(param.matrix, 1, median)
  if (ci == 0.05) {
    ci.lower <- apply(param.matrix, 1, function(x) { sort(x)[25] })
    ci.upper <- apply(param.matrix, 1, function(x) { sort(x)[975] })
  } else {
    if (ci == 0.01) {
      ci.lower <- apply(param.matrix, 1, function(x) { sort(x)[5] })
      ci.upper <- apply(param.matrix, 1, function(x) { sort(x)[995] })
    } else {
      if (ci == 0.1) {
        ci.lower <- apply(param.matrix, 1, function(x) { sort(x)[50] })
        ci.upper <- apply(param.matrix, 1, function(x) { sort(x)[950] })
      } else {
        stop("`ci` must be either 0.1, 0.05, or 0.01.")
      }
    }
  }

  ## Significance test
  is.sig <- logical(m)
  for (i in 1:m) {
    if (sign(ci.upper[i]) != sign(ci.lower[i])) {
      is.sig[i] <- FALSE
    } else {
      if (abs(brf.coef[i]) > abs((abs(ci.upper[i]) - abs(ci.lower[i]))/2)) {
        is.sig[i] <- TRUE
      } else {
        is.sig[i] <- FALSE
      }
    }
  }
    
  out <- data.frame(coef = brf.coef, significant = is.sig, ci.lower = ci.lower, ci.upper = ci.upper)
  rownames(out) <- colnames(p)
  if (sb) # close status bar (if TRUE)
    close(pb)
  attributes(out)$npar <- attributes(p)$npar
  attributes(out)$vnames <- vnames
  out
}

