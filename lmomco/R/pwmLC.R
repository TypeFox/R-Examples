"pwmLC" <-
function(x, threshold=NULL, nmom=5, sort=TRUE) {
  if(sort) x <- sort(x)

  if(is.null(threshold)) {
    warning("threshold is NULL")
    return(NULL)
  }
  T <- threshold
  x <- sapply(x,function(v) { if(v <= T) return(T); return(v)})

  n <- length(x)
  m <- length(x[x == T])
  observed.sample <- x[x > T]
  n.minus.m <- length(observed.sample)
  #  print(length(observed.sample))
  z <- pwm(observed.sample, nmom=nmom, sort=FALSE)
  #  print(z)
  Abetas <- z$betas

  Bbetas <- vector(mode="numeric",length=nmom)
  Bbetas <- rep(NA,length(Bbetas))

  #cat(c("DEBUG m.minus.1=",m.minus.1,"  n=",n,"\n"))

  for(r in seq(0,nmom-1)) {
    i <- r+1
    sumA <- 0
    sumB <- 0
    sumA <- sapply((m+1):n, function(j) { return(choose(j-1,r)*x[j]) })
    sumA <- sum(sumA)
    if(m > 0) { # avoid loop if no censored values
      sumB <- sapply(1:m, function(j) { return(choose(j-1,r)*T) } )
      sumB <- sum(sumB)
    }
    Bbetas[i] <- (sumA+sumB)/(n*choose(n-1,r))
  }

  zeta <- m/n

  z <- list(Aprimebetas=Abetas,
            Bprimebetas=Bbetas,
            source="pwmLC",
            threshold=T,
            zeta=zeta,
            numbelowthreshold=m,
            observedsize=(n-m),
            samplesize=n)
  return(z)

}
