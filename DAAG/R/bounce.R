"bounce" <-
function(y, d, log=FALSE)
{
  if(is.logical(log)&log)logbase <- 10
  else logbase <- log
  ord <- order(y)
  x <- y[ord]
  if(log)x <- log(x, base=logbase)
  n <- length(x)
  if(n <= 1)
    x
  else {
    i1 <- 1
    xnew <- x
    while(i1 < n) {
      x1 <- x[i1]
      i2 <- i1 + 1
      for(j in i2:n) {
        ## Identifies the longest contiguous sequence,
        ## starting at i1, of name slots.
        ## Because it is the longest such sequence, no
        ## subsequently sequence can overlap with it.
        nobounce <- TRUE
        jn <- n - j + i2
        dj <- x[jn] - x1
        dsought <- (jn - i1) * d
        if(dj < dsought) {
          jot <- (dsought - dj)/2
          for(k in i1:jn)
            xnew[k] <- x1 - jot + (k - i1) * d
          i1 <- jn + 1
          nobounce <- FALSE
          break
        }
      }
      if(nobounce)
        i1 <- i1 + 1
    }
    if(min(diff(xnew)) < d * 0.999) {
      n1 <- (1:(n - 1))[diff(xnew) < d]
      cat("Error in bounce().  Improperly separated points are:",
          fill = TRUE)
      cat(paste(n1, ":", n1 + 1, sep = ""), fill = TRUE)
      cat(paste(xnew[n1], ":", xnew[n1 + 1], sep = ""), fill
          = TRUE)
    }
    xnew
  }
  if(log)xnew <- logbase^xnew
  y[ord] <- xnew
  y
}

