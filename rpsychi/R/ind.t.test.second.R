ind.t.test.second <-
function(m, sd, n, unbiased=TRUE, correct=TRUE, sig.level=.05, digits=3){

m1 <- as.numeric(m[1]); m2 <- as.numeric(m[2])
sd1 <- as.numeric(sd[1]); sd2 <- as.numeric(sd[2])
n1 <- as.numeric(n[1]); n2 <- as.numeric(n[2])

  if(unbiased==FALSE){
     sd1 <- ssd2sd(n1,sd1)
     sd2 <- ssd2sd(n2,sd2)
  }


##(a) sample statistics
  samp.stat <- round(c(m1=m1, sd1=sd1, n1=n1, m2=m2, sd2=sd2, n2=n2), digits)

##(b) raw mean difference
  psi    <- m1 - m2
  dfw    <- n1 + n2 - 2
  sp.sq  <- ((n1-1) * sd1^2 + (n2-1) * sd2^2)/dfw 
  psi.std <- sqrt(sp.sq * (1/n1 + 1/n2))
  psi.lower <- psi + psi.std * qt(sig.level/2, dfw)
  psi.upper <- psi + psi.std * qt(sig.level/2, dfw, lower.tail=FALSE)  
  raw.difference <- round(c(mean.diff = psi, lower = psi.lower, upper = psi.upper, std=psi.std), digits)

##(c) standardized mean difference
  g     <- psi/sqrt(sp.sq)
  g.std <- sqrt(((g^2)/(2*dfw)) + (n1+n2)/(n1*n2))
  if(correct==TRUE){
    c.m <- (1 - 3/(4*dfw - 1))
    g <-  c.m * g
    g.std <- sqrt(((g^2)/(2*dfw)) + (n1+n2)/(n1*n2))
  }  
  g.lower <- g + g.std * qnorm(sig.level/2)
  g.upper <- g + g.std * qnorm(sig.level/2, lower.tail=FALSE)
  standardized.difference <- round(c(es=g, lower=g.lower, upper=g.upper, std=g.std), digits)
  


##(d) Power Analysis
  c.delta <- c(.2, .5, .8)
  criterion.power <- round(power.t(sig.level=sig.level,delta=c.delta,n1=n1,n2=n2), digits)
  names(criterion.power) <- c("small", "medium", "large")

## output
  output <- list(samp.stat=samp.stat, raw.difference=raw.difference, standardized.difference=standardized.difference, power=criterion.power)
  return(output)
}

