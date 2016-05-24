



# Initial values (if not provided)


inits<-function(y,x,init,q,m) {
p <- ncol(x)

fit <- lm( y ~ -1 + x )
if (missing(init) || is.na(match("beta", names(init) ) ) )
{
  beta <- as.vector(fit$coef);   names(beta) <- NULL
}
else {
  beta <- as.vector(init$beta); names(beta) = NULL
}

if (missing(init) || is.na(match("sigma", names(init) ) ) )
  sigma2 <- (summary(fit)$sigma)^2
else
  sigma2 <- init$sigma^2

if (missing(init) || is.na(match("Delta", names(init) ) ) )
  Delta <- diag(rep(1,q))


else
  Delta = init$Delta

if (missing(init) || is.na(match("bi", names(init) ) ) )
  b <- matrix(rep(0,m*q), ncol=m)
else
  b = init$bi
 
 return(list(beta=beta,sigma2=sigma2,Delta=Delta,b=b)) 
 
  }
  
  