"formtomat" <-
  function(a, xvar) {

  d <- strsplit(a, "\\+")[[1]]
  y <- as.list(rep(NA, times = (length(d) + 1)))
  for (i in 1:length(xvar)) {
    w <- regexpr(xvar[i], d)
    for (j in 1:length(w)) {
      if (w[j] != -1) {
        if (is.na(y[[j+1]])) {
          y[[j+1]] <- i
        }
        else {
          y[[j+1]] <- c(y[[j+1]], i)
        }
      }
    }
    w2 <- regexpr("\\^2", d)
    for (j in 1:length(w2)) {
      if ((w2[j] != -1) & (w[j] != -1)) {
        y[[j+1]] <- c(y[[j+1]], i)
      }
    }
  }

  yp <- as.list(rep(NA, times = length(xvar)))
  for (i in 1:length(xvar)) {
    yp[[i]] <- as.list(rep(NA, times = (length(d)+1)))
  }
  np <- as.list(rep(NA, times = length(xvar)))
  for (i in 1:length(xvar)) {
    np[[i]] <- rep(1, times = (length(d)+1))
  }
  for (j in 1:length(xvar)) {
    np[[j]][1] <- 0
  }
  for (i in 1:length(d)) {
    a.p <- parse(text = d[i])
    for (j in 1:length(xvar)) {
      str <- deparse(D(a.p, xvar[j]))
      w <- regexpr("[0-8]", str)
      if (w != -1) {
        np[[j]][i+1] <- as.numeric(substring(str, w, w+1))
      }
      for (k in 1:length(xvar)) {
        w <- regexpr(xvar[k], str)
        if (w != -1) {
          if (is.na(yp[[j]][[i+1]])) {
            yp[[j]][[i+1]] <- k
          }
          else {
            yp[[j]][[i+1]] <- c(y[[j]][[i+1]], k)
          }
        }
      }
    }
  }

  return(list(ind = y, ind.n = np, ind.p = yp))
  
}
