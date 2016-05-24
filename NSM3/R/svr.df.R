svr.df<-function (z, delta,lambda.hat=0.001, alpha = 3, npoints=2053) 
{
  X <- cbind(delta,z)
  lambda.j <- 1
  n <- dim(X)[1]
  
  r = sum(delta==1)
  # Reorder
  X <- X[rev(order(X[, 1])), ]
  Z.star <- sort(X[X[, 1] == 0, 2])
  Z <- sort(X[, 2])    
  m = length(Z.star) 
  
  
  product <- 1
  # Get factors for product
  for(j in 1:m)
  {
    N.n <- length(Z[Z >= Z.star[j]])
    numer <- alpha * pexp(Z.star[j], rate = lambda.hat, lower.tail=F) + N.n
    denom <- numer - lambda.j
    product <- c(product, (numer / denom))
  }
  
  x <- seq(0, max(Z), length=npoints)
  
  F.hat <- numeric(0)
  for(l in 1:(m + 1))  # "ell"
  {
    # Segment x with Z.star
    if(l == 1) x.start <- 1
    if(l > 1) x.start <- x.start.next
    if(l <= m) x.end <- length(x[x < Z.star[l]])
    if(l == m + 1) x.end <- length(x)
    x.start.next <- x.end + 1
    x.l <- x[x.start:x.end]
    
    # Evaluate at x.l
    N.n <- numeric(0)
    for(i in x.l) N.n <- c(N.n, length(Z[Z > i]))
    F.hat.l <- alpha * pexp(x.l, rate=lambda.hat, lower.tail=F) + N.n
    F.hat.l <- F.hat.l / (alpha + n)
    F.hat.l <- F.hat.l * prod(product[1:l])
    F.hat <- c(F.hat, F.hat.l)
  }
  
  plot(x, F.hat, type="l")
  #plot(survfit(Surv(hodgkins.affected$days, hodgkins.affected$relapse)~1, conf.type="none"), 
  #     lty=2, xlab="Survival Time in Days (t)",
  #     ylab="Estimated Probability of No Relapse to Time t")  
  #print(F.hat)    
  for(l in 1:26)
  {
    # Segment x with Z
    if(l == 1) x.start <- 1
    if(l > 1) x.start <- x.start.next
    if(l <= 25) x.end <- length(x[x < Z[l]])
    if(l == 26) x.end <- length(x)
    x.start.next <- x.end + 1
    if(l == 1) lines(x[x.start:x.end], F.hat[x.start:x.end],  
                     xlim=c(0, max(Z)), ylim=c(0, max(F.hat)), type="l",
                     xlab="Survival Time in Days (t)",
                     ylab="Estimated Probability of No Relapse to Time t")
    if(i > 1) lines(x[x.start:x.end], F.hat[x.start:x.end])
  }
  legend("topright", lty=c(1:2), legend=c("Susarla-Van Ryzin", "Kaplan-Meier"))
  
  
  
  list(x, F.hat)
}
