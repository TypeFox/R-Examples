log2 <- function(x) ifelse(x > 0,log(x), 0)

mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))

# compute the bivariate normal cdf and its derivates
mypbivnorm <- function(x1, x2, rho = 0){
  if (is.null(x1) &&  is.null(x2)) result <- list(f = 1, a = 0, b = 0, rho = 0, rho.rho = 0)
  if (is.null(x1) && !is.null(x2)) result <- list(f = pnorm(x2), a = 0, b = dnorm(x2), rho = 0, rho.rho = 0)
  if (is.null(x2) && !is.null(x1)) result <- list(f = pnorm(x1), a = dnorm(x1), b = 0, rho = 0, rho.rho = 0)
  if (!is.null(x1) && !is.null(x2)){
    f <- pbivnorm(x1, x2, rho)
    b <- dnorm(x2) * pnorm( (x1 - rho * x2) / sqrt(1 - rho^2) )
    a <- dnorm(x1) * pnorm( (x2 - rho * x1) / sqrt(1 - rho^2) )
    eps <- 1E-07
    rho <- (pbivnorm(x1, x2, rho + 1E-07) - pbivnorm(x1, x2, rho)) / 1E-07
    result <- list(f = f, a = a, b = b, rho = rho)
  }
  result
}

# a function to construct block-diagonal matrix (can't remember from
# which package it is borrowed)
bdiag <- function(...){
  if (nargs() == 1)
    x <- as.list(...)
  else
    x <- list(...)
  n <- length(x)
  if(n == 0) return(NULL)
  x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
              stop("Zero-length component in x"))
  d <- array(unlist(lapply(x, dim)), c(2, n))
  rr <- d[1,]
  cc <- d[2,]
  rsum <- sum(rr)
  csum <- sum(cc)
  out <- array(0, c(rsum, csum))
  ind <- array(0, c(4, n))
  rcum <- cumsum(rr)
  ccum <- cumsum(cc)
  ind[1,-1] <- rcum[-n]
  ind[2,] <- rcum
  ind[3,-1] <- ccum[-n]
  ind[4,] <- ccum
  imat <- array(1:(rsum * csum), c(rsum, csum))
  iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
                                               (y[3]+1):y[4]], imat=imat)
  iuse <- as.vector(unlist(iuse))
  out[iuse] <- unlist(x)
  return(out)
} 

# print a summary of the non-linear opimization
print.est.stat <- function(x, ...){
  et <- x$elaps.time[3]
  i <- x$nb.iter[1]
  halton <- x$halton
  method <- x$method
  if (!is.null(x$type) && x$type != "simple"){
    R <- x$nb.draws
    cat(paste("Simulated maximum likelihood with", R, "draws\n"))
  }
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  cat(paste(method, "method\n"))
  tstr <- paste(h, "h:", m, "m:", s, "s", sep="")
  cat(paste(i,"iterations,",tstr,"\n"))
  if (!is.null(halton)) cat("Halton's sequences used\n")
  if (!is.null(x$eps)) cat(paste("g'(-H)^-1g =", sprintf("%5.3G", as.numeric(x$eps)),"\n"))
  if (is.numeric(x$code)){
    msg <- switch(x$code,
                  "1" = "gradient close to zero",
                  "2" = "successive fonction values within tolerance limits",
                  "3" = "last step couldn't find higher value",
                  "4" = "iteration limit exceeded"
                  )
    cat(paste(msg, "\n"))
  }
  else cat(paste(x$code, "\n"))
}


