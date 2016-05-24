"binom.exact" <-
  function(x, n, conf.level=.95) {
    xnc <- cbind(x,n,conf.level)
    lower <- numeric(nrow(xnc))
    upper <- numeric(nrow(xnc))
    for(i in 1:nrow(xnc)){
      ci <- binom.test(x=xnc[i,1], n=xnc[i,2],
                       conf.level=xnc[i,3])$conf.int
      lower[i] <- ci[1]
      upper[i] <- ci[2]
    }
    data.frame(x = x,
               n = n,
               proportion = x/n,
               lower = lower,
               upper = upper,
               conf.level = conf.level
               )
}
