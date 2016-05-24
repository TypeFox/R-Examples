killworth <- function(dat, known, N) {
  n <- dim(dat)[1]
  indices <- 1:length(known)
  indices.k <- (length(known)+1):(dim(dat)[2])
  
  NSUM.d.back <- vector(length=n)
  for(j in 1:n)
  {
    NSUM.d.back[j] <- sum(dat[j, indices])*(N/sum(known))
  }
  return(N*(apply(as.matrix(dat[,indices.k]),2,sum)/sum(NSUM.d.back)))
}

killworth.start <- function(dat, known, N) {
  n <- dim(dat)[1]
  indices <- 1:length(known)
  indices.k <- (length(known)+1):(dim(dat)[2])
  
  NSUM.d.back <- vector(length=n)
  for(j in 1:n)
  {
    NSUM.d.back[j] <- sum(dat[j, indices])*(N/sum(known))
  }
  
  NK.start <- N*(apply(as.matrix(dat[,indices.k]),2,sum)/sum(NSUM.d.back))
  
  d.start <- NSUM.d.back
  for(j in 1:n)
  {
    if(d.start[j] < max(dat[j,]))
    {
      d.start[j] <- max(dat[j,]) + 1
    }
  }	
  
  if(sort(NSUM.d.back)[1] == 0)
  {
    zero.index <- which(NSUM.d.back == 0)
    d.start[zero.index] <- 1
    NSUM.d.back <- NSUM.d.back[-zero.index]
  }
  lnorm.start <- as.vector(fitdistr(NSUM.d.back, "log-normal")$estimate)
  mu.start <- lnorm.start[1]
  sigma.start <- lnorm.start[2]
  
  return(list(NK.start=NK.start, d.start=d.start, mu.start=mu.start, sigma.start=sigma.start))
}
