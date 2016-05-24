"fnfam" <-
function(n,power,alpha,p1,h)
 {
  nhp <- n*2*h
  s <- sqrt(p1*(1-p1)/nhp)
 pb <- p1 - s*qnorm(power)
 f <- 2*n*2*h*(pb*log(pb) + (1-pb)*log(1-pb) - log(0.5)) - (qnorm(alpha))^2 
 return(f)
 }

