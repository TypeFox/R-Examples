# Linearizable C-E1E1
# Common variance for both segments

con.search.E1E1.C <- function(x, y, n, jlo, jhi)
{
  fjk <- matrix(0, n)
  fxy <- matrix(0, jhi - jlo + 1)
  
  jkgrid <- expand.grid(jlo:jhi)
  res <- data.frame(j = jkgrid,
                    k.ll = apply(jkgrid, 1, con.parmsFUN.E1E1.C, x = x, 
                                 y = y, n = n))
  
  fxy <- matrix(res$k.ll, nrow = jhi-jlo+1)
  rownames(fxy) <- jlo:jhi
  
  z <- findmax(fxy)
  jcrit <- z$imax + jlo - 1
  list(jhat = jcrit, value = max(fxy))
}

con.parmsFUN.E1E1.C <- function(j, x, y, n){
  a <- con.parms.E1E1.C(x,y,n,j,1)
  nr <- nrow(a$theta)
  est <- a$theta[nr,  ]
  b<-con.est.E1E1.C(x[j],est)
  s2<-1/b$eta1
  return(p.ll.C(n, j, s2))
}

con.parms.E1E1.C <- function(x,y,n,j0,e10){
  
  th <- matrix(0,100,4)
  
  # Iteration 0
  
  th[1,1] <- e10
  bc <- beta.calc.E1E1.C(x,y,n,j0,e10)
  th[1,2:4] <- bc$B
  
  # Iterate to convergence (100 Iter max)
  
  for (iter in 2:100){
    m <- iter-1
    ec <- eta.calc.E1E1.C(x,y,n,j0,th[m,2:4])
    th[iter,1] <- ec$eta1
    bc <- beta.calc.E1E1.C(x,y,n,j0,ec$eta1)
    th[iter,2:4] <- bc$B
    theta <- th[1:iter,]
    #delta <- abs(th[iter,]-th[m,])
    delta <- abs(th[iter,]-th[m,])/th[m,]
    if( (delta[1]<.001) & (delta[2]<.001) & (delta[3]<.001)
        & (delta[4]<.001) )
      break
  }
  list(theta=theta)
}
con.est.E1E1.C <- function(xj, est)
{
  eta1 <- est[1]
  a0 <- est[2]
  a1 <- est[3]
  b1 <- est[4]
  b0 <- a0 + (a1 - b1) * exp(xj)
  list(eta1 = eta1, a0 = a0, a1 = a1, b0 = b0, b1 = b1)
}
con.vals.E1E1.C <- function(x, y, n, j)
{
  a <- con.parms.E1E1.C(x, y, n, j, 1)
  nr <- nrow(a$theta)
  est <- a$theta[nr,  ]
  b <- con.est.E1E1.C(x[j], est)
  eta <- c(b$eta1)
  beta <- c(b$a0, b$a1, b$b0, b$b1)
  tau <- x[j]
  list(eta = eta, beta = beta, tau = tau)
}

p.ll.C <- function(n, j, s2){
  q1 <- n * log(sqrt(2 * pi))
  q2 <- 0.5 * (n) * (1 + log(s2))
  - (q1 + q2)
}

findmax <-function(a){
  maxa<-max(a)
  imax<- which(a==max(a),arr.ind=TRUE)[1]
  jmax<-which(a==max(a),arr.ind=TRUE)[2]
  list(imax = imax, jmax = jmax, value = maxa)
}

beta.calc.E1E1.C <- function(x, y, n, j, e1)
{
  aa <- wmat.E1E1.C(x, y, n, j, e1)
  W <- aa$w
  bb <- rvec.E1E1.C(x, y, n, j, e1)
  R <- bb$r
  beta <- solve(W, R)
  list(B = beta)
} 
eta.calc.E1E1.C <- function(x, y, n, j, theta)
{
  jp1 <- j + 1
  a0 <- theta[1]
  a1 <- theta[2]
  b1 <- theta[3]
  b0 <- a0 + (a1 - b1) * exp(x[j])
  rss1 <- sum((y[1:j] - a0 - a1 * exp(x[1:j]))^2)
  rss2 <- sum((y[jp1:n] - b0 - b1 * exp(x[jp1:n]))^2)
  e1 <- (n)/(rss1 + rss2)
  list(eta1 = e1)
}
wmat.E1E1.C <- function(x, y, n, j, e1)
{
  W <- matrix(0, 3, 3)
  jp1 <- j + 1
  W[1, 1] <- e1 * n 
  W[1, 2] <- e1 * (sum(exp(x[1:j])) + (n - j) * exp(x[j]))
  W[1, 3] <- e1 * sum( exp(x[jp1:n]) - exp(x[j]))
  
  W[2, 2] <- e1 * sum(exp(x[1:j])^2) + e1 * (n-j) * exp(x[j])^2
  W[2, 3] <- e1 * exp(x[j]) * sum( exp(x[jp1:n]) - exp(x[j]) )
  
  W[3, 3] <-e1 * sum(  ( exp(x[jp1:n]) - exp(x[j]) )^2 )
  
  W[2, 1] <- W[1, 2]
  W[3, 1] <- W[1, 3]
  W[3, 2] <- W[2, 3]
  
  list(w = W)
}
rvec.E1E1.C <- function(x, y, n, j, e1)
{
  R <- array(0, 3)
  jp1 <- j + 1
  y1j <- sum(y[1:j])
  yjn <- sum(y[jp1:n])
  xy1j <- sum(exp(x[1:j]) * y[1:j])
  xyjn <- sum(exp(x[jp1:n]) * y[jp1:n])
  
  R[1] <- e1 * (y1j + yjn)
  R[2] <- e1 * (xy1j + exp(x[j]) * yjn)
  R[3] <- e1 * (xyjn - yjn * exp(x[j]))
  list(r = R)
}


