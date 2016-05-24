lpjk <-
  function(jk, x, y, n, sigma2, tau2)
  {  
    j = jk[1]
    k = jk[2]
    g <- rss(x, y, n, j, k)
    rss1 <- g$S1
    rss2 <- g$S2
    rss3 <- g$S3
    x1 <- ((k - j)/2) * log(sigma2/tau2)
    x2 <- (rss1 + rss3)/(2 * sigma2)
    x3 <- rss2/(2 * tau2)
    x1 - x2 - x3
  }
