# Linear C-LQ
# Common variance for both segments

con.search.LQ.C <- function(x, y, n, jlo, jhi)
{
  fjk <- matrix(0, n)
  fxy <- matrix(0, jhi - jlo + 1)
  
  jkgrid <- expand.grid(jlo:jhi)
  res <- data.frame(j = jkgrid,
                    k.ll = apply(jkgrid, 1, con.parmsFUN.LQ.C, x = x, 
                                 y = y, n = n))
  
  fxy <- matrix(res$k.ll, nrow = jhi-jlo+1)
  rownames(fxy) <- jlo:jhi
  
  z <- findmax(fxy)
  jcrit <- z$imax + jlo - 1
  list(jhat = jcrit, value = max(fxy))
}

con.parmsFUN.LQ.C <- function(j, x, y, n){
  a <- con.parms.LQ.C(x,y,n,j,1)
  nr <- nrow(a$theta)
  est <- a$theta[nr,  ]
  b<-con.est.LQ.C(x[j],est)
  s2<-1/b$eta1
  return(p.ll.C(n, j, s2))
}

con.parms.LQ.C <- function(x,y,n,j0,e10){
  
  th <- matrix(0,100,5)
  
  # Iteration 0
  
  th[1,1] <- e10
  bc <- beta.calc.LQ.C(x,y,n,j0,e10)
  th[1,2:5] <- bc$B
  
  # Iterate to convergence (100 Iter max)
  
  for (iter in 2:100){
    m <- iter-1
    ec <- eta.calc.LQ.C(x,y,n,j0,th[m,2:5])
    th[iter,1] <- ec$eta1
    bc <- beta.calc.LQ.C(x,y,n,j0,ec$eta1)
    th[iter,2:5] <- bc$B
    theta <- th[1:iter,]
    #delta <- abs(th[iter,]-th[m,])
    delta <- abs(th[iter,]-th[m,])/th[m,]
    if( (delta[1]<.001) & (delta[2]<.001) & (delta[3]<.001)
        & (delta[4]<.001) & (delta[5]<.001) )
      break
  }
  list(theta=theta)
}
con.est.LQ.C <- function(xj, est)
{
  eta1 <- est[1]
  a0 <- est[2]
  a1 <- est[3]
  b1 <- est[4]
  b2 <- est[5]
  b0 <- a0 + (a1 - b1) * xj - b2 * xj^2
  list(eta1 = eta1, a0 = a0, a1 = a1, b0 = b0, b1 = b1, b2 = b2)
}

con.vals.LQ.C <- function(x, y, n, j)
{
  a <- con.parms.LQ.C(x, y, n, j, 1)
  nr <- nrow(a$theta)
  est <- a$theta[nr,  ]
  b <- con.est.LQ.C(x[j], est)
  eta <- c(b$eta1)
  beta <- c(b$a0, b$a1, b$b0, b$b1, b$b2)
  tau <- x[j]
  list(eta = eta, beta = beta, tau = tau)
}

p.ll.C <- function(n, j, s2){
  q1 <- n * log(sqrt(2 * pi))
  q2 <- 0.5 * n * (1 + log(s2))
  - (q1 + q2)
}

findmax <-function(a){
  maxa<-max(a)
  imax<- which(a==max(a),arr.ind=TRUE)[1]
  jmax<-which(a==max(a),arr.ind=TRUE)[2]
  list(imax = imax, jmax = jmax, value = maxa)
}

beta.calc.LQ.C <- function(x, y, n, j, e1)
{
  aa <- wmat.LQ.C(x, y, n, j, e1)
  W <- aa$w
  bb <- rvec.LQ.C(x, y, n, j, e1)
  R <- bb$r
  beta <- solve(W, R)
  list(B = beta)
} 

eta.calc.LQ.C <- function(x, y, n, j, theta)
{
  jp1 <- j + 1
  a0 <- theta[1]
  a1 <- theta[2]
  b1 <- theta[3]
  b2 <- theta[4]
  b0 <- a0 + (a1 - b1) * x[j] - b2 * x[j]^2
  rss1 <- sum((y[1:j] - a0 - a1 * x[1:j])^2)
  rss2 <- sum((y[jp1:n] - b0 - b1 * x[jp1:n] - b2 * x[jp1:n]^2)^2)
  e1 <- n/(rss1 + rss2)
  list(eta1 = e1)
}

wmat.LQ.C <- function(x, y, n, j, e1)
{
  W <- matrix(0, 4, 4)
  jp1 <- j + 1
  W[1, 1] <- e1 * n 
  W[1, 2] <- e1 * sum(x[1:j]) + e1 * (n - j) * x[j]
  W[1, 3] <- e1 * sum(x[jp1:n] - x[j]) 
  W[1, 4] <- e1 * sum(x[jp1:n]^2 - x[j]^2)
  
  W[2, 2] <- e1 * sum(x[1:j]^2) + e1 * (n - j) * x[j]^2
  W[2, 3] <- e1 * x[j] * sum(x[jp1:n] - x[j])
  W[2, 4] <- e1 * x[j] * sum(x[jp1:n]^2 - x[j]^2)
  
  W[3, 3] <- e1 * sum( (x[jp1:n] - x[j])^2)
  W[3, 4] <- e1 * sum( (x[jp1:n]^2 - x[j]^2) * (x[jp1:n] - x[j]) )
  
  W[4, 4] <- e1 * sum((x[jp1:n]^2 - x[j]^2)^2)
    
  W[2, 1] <- W[1, 2]
  W[3, 1] <- W[1, 3]
  W[4, 1] <- W[1, 4]
  W[3, 2] <- W[2, 3]
  W[4, 2] <- W[2, 4]
  W[4, 3] <- W[3, 4]
  
  list(w = W)
}

rvec.LQ.C <- function(x, y, n, j, e1)
{
  R <- array(0, 4)
  jp1 <- j + 1
  y1j <- sum(y[1:j])
  yjn <- sum(y[jp1:n])
  xy1j <- sum(x[1:j] * y[1:j])
  xyjn <- sum(x[jp1:n] * y[jp1:n])
  x2yjn <- sum(x[jp1:n]^2 * y[jp1:n])
  
  R[1] <- e1 * (y1j + yjn)
  R[2] <- e1 * (xy1j + x[j] * yjn)
  R[3] <- e1 * (xyjn - yjn * x[j])
  R[4] <- e1 * (x2yjn - x[j]^2 * yjn)
  list(r = R)
}


