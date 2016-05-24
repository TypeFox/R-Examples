dep.t.test.second <- function(m, sd, n, corr, unbiased=TRUE, sig.level=.05, digits=3){
  m1 <- as.numeric(m[1]); m2 <- as.numeric(m[2])
  sd1 <- as.numeric(sd[1]); sd2 <- as.numeric(sd[2])

  if(unbiased==FALSE){
     sd1 <- ssd2sd(n, sd1)
     sd2 <- ssd2sd(n, sd2)
  }
  

##(a) sample statistics
  samp.stat <- round(c(m1=m1, sd1=sd1, m2=m2, sd2=sd2, n=n, corr=corr), digits)


##(b) raw mean difference
  psi    <- m1 - m2
  dfw    <- n - 1
  cov <- corr * sd1 * sd2 
  psi.std <- sqrt((sd1^2 + sd2^2 - 2 * cov) /n)
  psi.lower <- psi + psi.std * qt(sig.level/2, dfw)
  psi.upper <- psi + psi.std * qt(sig.level/2, dfw, lower.tail=FALSE)
  raw.difference <- round(c(mean.diff = psi, lower = psi.lower, upper = psi.upper, std=psi.std), digits)
  

##(c) standardized mean difference
  sp <- sqrt((dfw * sd1^2 + dfw * sd2^2)/ (dfw*2))  #4.5
  g  <- psi / sp
 	g.std <- sqrt(g^2 /(2 * (n-1)) + 2 * (1 - corr) / n)
 	g.lower <-  g + g.std * qnorm(sig.level/2,lower.tail=TRUE)
 	g.upper <-  g + g.std * qnorm(sig.level/2,lower.tail=FALSE)
  standardized.difference <- round(c(es=g, lower=g.lower, upper=g.upper, std=g.std), digits)  
  
  
  output <- list(samp.stat=samp.stat, raw.difference=raw.difference, standardized.difference=standardized.difference)
  return(output) 
}
