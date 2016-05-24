## Generate samples from mle log-concave density
'rlcd' <- function (n = 1, lcd, method = c("Independent","MH")) {
  method <- match.arg(method)
  triang <- lcd$triang
  x <- lcd$x
  logMLE <- lcd$logMLE
  nrows <- nrow(triang)
  d <- ncol(x)
  prob <- rep(0, nrows)
  for (j in 1:nrows) {
    prob[j] <- lcd$detA[j] * JAD(lcd$logMLE[triang[j, ]])
  }
  prob <- prob/sum(prob)
  samples <- matrix(0, nrow = n, ncol = d)
  
  ## pick a simplex for each sample
  simp <- sample(1:nrows, n, prob = prob, replace = TRUE)
  if (method == "MH"){
    px = 0
    qx = 0
    while (sum(samples[1, ] == 0)) {
      ## generate point on unit simplex
      w <- rexp(d + 1)
      w <- w/sum(w)
      y <- logMLE[triang[simp[1], ]]
      ##evaluate at the corresponding point
      fw <- exp(y %*% w)
      maxfx <- max(exp(y))
      u <- runif(1)
      if (u < fw/maxfx) {
        samples[1, ] <- w %*% x[triang[simp[1], ], ]
        px = fw / prob[simp[1]] * lcd$detA[simp[1]]
      }
    }
    for (i in 2:n) {
      w <- rexp(d + 1)
      w <- w/sum(w)
      y <- logMLE[triang[simp[i], ]]
      fw <- exp(y %*% w)
      qx = fw / prob[simp[i]] * lcd$detA[simp[i]]
      u <- runif(1)
      if (u < min(qx/px,1)) {
        samples[i, ] <- w %*% x[triang[simp[i], ], ]
        px = qx
      }
      else samples[i, ] = samples[i-1, ]       
    }
  }
  else {
    for (i in 1:n) {
      while (sum(samples[i, ] == 0)) {
        w <- rexp(d + 1)
        w <- w/sum(w)
        y <- logMLE[triang[simp[i], ]]
        fw <- exp(y %*% w)
        maxfx <- max(exp(y))
        u <- runif(1)
        if (u < fw/maxfx) samples[i, ] <- w %*% x[triang[simp[i], ], ]
      }
    }
  }    
  return(samples)
}














