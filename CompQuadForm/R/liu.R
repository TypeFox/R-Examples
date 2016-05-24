liu <- function(q, lambda, h = rep(1,length(lambda)), delta = rep(0,length(lambda))) {

  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
 
  c1 <- sum(lambda*h) + sum(lambda*delta)

  c2 <- sum(lambda^2*h) + 2*sum(lambda^2*delta)

  c3 <- sum(lambda^3*h) + 3*sum(lambda^3*delta)

  c4 <- sum(lambda^4*h) + 4*sum(lambda^4*delta)
  
  s1 <- c3/(c2^(3/2))

  s2 <- c4/c2^2

  muQ <- c1

  sigmaQ <- sqrt(2*c2)

  tstar <- (q-muQ)/sigmaQ

  if (s1^2>s2) {

    a <- 1/(s1-sqrt(s1^2-s2))

    delta <- s1*a^3-a^2

    l <- a^2-2*delta

  } else {

    a <- 1/s1
    
    delta <- 0

    l <- c2^3/c3^2

  }

  muX <- l+delta

  sigmaX <- sqrt(2)*a
  
  Qq <- pchisq(tstar*sigmaX+muX,df=l,ncp=delta,lower.tail=FALSE)

  return(Qq)

}
