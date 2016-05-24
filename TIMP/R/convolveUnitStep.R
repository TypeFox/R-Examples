convolveUnitStep <- function(k, x, irfpar) {

  m <- matrix(0, nrow=length(x), ncol=length(k))
  
  # save a little time by storing the unitstep results
  step1 <- matrix(0, nrow=length(x), ncol=length(irfpar))
  for(j in 1:length(irfpar)) 
    step1[,j] <- unitstep(x-irfpar[j])
  for(i in 1:length(k)) {
    dec <- exp(k[i]*x)
    for(j in 1:length(irfpar)) {
      s1 <- if(j%%2==0) -1 else 1
      s2 <- -1 * s1
      decj <- exp(k[i]*irfpar[j])
      
      m[,i] <- m[,i] +
        (( s1 * dec + s2 * decj) * step1[,j]) +
          ((-1 + decj) * unitstep(- irfpar[j]))     
      
    }
    m[,i] <- m[,i] * exp(-k[i]*x) * (1/k[i]) 
  }
  m
}
