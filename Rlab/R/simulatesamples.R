"simulatesamples" <-function(nrep, n, DIST, data = FALSE, ...){
  cat("Simulated", nrep, " data sets each with ", n, "obserations", fill = TRUE)
  temp <- matrix(DIST(n * nrep, ...), ncol = n, nrow = nrep)
  xbar <- c(temp %*% rep(1./n, n))
  s <- c(sqrt((1./(n - 1.)) * ((temp - matrix(xbar, nrow = nrep, ncol = n))^2.) %*% rep(1., n)))
  if((n * nrep > 1000.) && data) {
    cat("don't forget that the matrix of simluated samples is a fairly large dataset", fill = TRUE)
    cat("( ", n * nrep, " values total) and you should delete this data set when you are", fill = TRUE)
    cat("finished analyzing it.", fill = TRUE)
    }
  if(data) return(list(dat = temp, m = xbar, s = s))
  else return(list(m = xbar, s = s))
  }