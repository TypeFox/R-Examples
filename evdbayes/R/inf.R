## COMPUTE POSTERIOR QUANTILES IN UPPER TAIL
"mc.quant" <-
  function(post, p, lh = c("gev", "gpd"))
{
  nc <- length(p)
  nr <- nrow(post)
  dn <- list(rownames(post), p)
  mat <- matrix(0, ncol = nc, nrow = nr, dimnames = dn)
  loc <- post[,"mu"]
  scale <- post[,"sigma"]
  shape <- post[,"xi"]
  
  if (lh == "gev")
    y <- -log(p)
  else
    y <- 1-p

  for(i in 1:nc)
    mat[,i] <- ifelse(shape,
                      loc + scale * (y[i]^(-shape) - 1)/shape,
                      loc - scale * log(y[i]))
  
  drop(mat)
}

## POSTERIOR RETURN LEVEL PLOT
"rl.pst" <-
  function(post, npy, lh = c("gev", "gpd"), ci = 0.9, lty = c(2,1),
           col = c(2,1), xlab = "return period",
           ylab = "return level",  ...)
  {
    if (missing(npy) && lh == "gpd")
      stop("``npy'' should be present with a ``gpd'' likelihood")

    if (lh == "gev")
      npy <- 1
    
    rps <- c(1/npy + 0.001, 10^(seq(0,4,len=20))[-1])
    p.upper <- 1 - 1/(npy * rps)
    mat <- mc.quant(post = post, p = p.upper, lh = lh) 
    mat <- t(apply(mat, 2, quantile, probs = c((1-ci)/2, 0.5, (1+ci)/2)))
    matplot(rps, mat, log = "x", type = "l",
            xlab = xlab, ylab = ylab, lty = lty, col = col, ...)
    invisible(list(x = rps, y = mat))
  }

"rl.pred" <-
  function(post, qlim, npy, lh = c("gev", "gpd"), period = 1, lty = 1,
           col = 1, xlab = "return period",
           ylab = "return level", ...)
{

  if (missing(npy) && lh == "gpd")
    stop("``npy'' should be present with a ``gpd'' likelihood")
  
  if (lh == "gev")
    npy <- 1
  
  np <- length(period)
  p.upper <- matrix(0, nrow = 25, ncol = np)
  qnt <- seq(qlim[1], qlim[2], length = 25)
  
  for(i in 1:25) {
    p <- (qnt[i] - post[,"mu"])/post[,"sigma"]
    p <- ifelse(post[,"xi"],
                exp( - pmax((1 + post[,"xi"] * p),0)^(-1/post[,"xi"])),
                exp(-exp(-p)))
    for(j in 1:np)
      p.upper[i,j] <- 1-mean(p^period[j])
  }

  if (lh == "gpd")
    p <- 1 + log(p)
  
  if(any(p.upper == 1))
    stop("lower q-limit is too small")
  if(any(p.upper == 0))
    stop("upper q-limit is too large")
  matplot(1/(npy * p.upper ), qnt, log = "x", type = "l", lty = lty,
          col = col, xlab = xlab, ylab = ylab, ...)
  invisible(list(x = 1/(npy * p.upper ), y = qnt))
}





