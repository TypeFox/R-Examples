distIneqMassartPlot <- function(densFn = "norm", param = NULL,
                                nSamp = 50, n = 100, ...){
  ## Using repeated sampling to get the LHS probability
  CALL <- match.call()
  dfun <- match.fun(paste("d", densFn, sep = ""))
  pfun <- match.fun(paste("p", densFn, sep = ""))
  rfun <- match.fun(paste("r", densFn, sep = ""))
  ## RHS bound
  RHS <- seq(0, 2, length.out = 200)

  ## Value of t for that probability bound
  tVal <- sqrt(log(RHS/2)/(-2*n))

  ## Number of RHS <= LHS
  m <- 0
  pm <- 0

  for (j in 1:200){
    for (i in 1:nSamp){
      if(is.null(param))
        x <- rfun(n = n, ...)
      else x <- rfun(n = n, param = param)
      ## calculate sup of empirical from true distribution function
      tx <- table(x)
      xi <- as.numeric(names(tx))
      if(is.null(param))
        f <- pfun(xi, ...)
      else f <- pfun(xi, param = param)
      fhat <- cumsum(tx)/n
        diff <- max(abs(fhat-f))
      pm <- ifelse(diff > tVal[j], pm + 1, pm)
    }
    m[j] = pm
    pm = 0
  }
  prob = m/nSamp

  ## Smooth the curve
  smooth <- 0
  probs <- 0

  for (k in 1:191){
    smooth[k] = mean(prob[k + 0:9])
  }

  for (m in 1:9){
    probs[m] = mean(prob[1:m])
  }

  smooth = c(probs,smooth)

  ## Plot the probability and its bound
  plot(tVal, RHS, type = "l",
       main = paste("The Massart Inequality for p", densFn, sep = ""),
       xlab = "t", ylab = "Probability and bound")
  lines(tVal, smooth, col = "red")
  invisible(NULL)
}
