# C-E2E2
# Different variances for two segments

con.search.E2E2.D <- function(x, y, n, jlo, jhi)
{
  fjk <- matrix(0, n)
  fxy <- matrix(0, jhi - jlo + 1)
  
  jkgrid <- expand.grid(jlo:jhi)
  res <- data.frame(j = jkgrid,
                    k.ll = apply(jkgrid, 1, con.parmsFUN.E2E2.D, x = x, 
                                 y = y, n = n))
  
  fxy <- matrix(res$k.ll, nrow = jhi-jlo+1)
  rownames(fxy) <- jlo:jhi
  
  z <- findmax(fxy)
  jcrit <- z$imax + jlo - 1
  list(jhat = jcrit, value = max(fxy))
}

con.parmsFUN.E2E2.D <- function(j, x, y, n){
  a <- con.parms.D.E2E2(x,y,n,j,1,1)
  nr <- nrow(a$theta)
  est <- a$theta[nr,  ]
  b<-con.est.D.E2E2(x[j],est)
  s2<-1/b$eta1
  t2<-1/b$eta2
  return(p.ll.D(n, j, s2, t2))
}

con.parms.D.E2E2 <-
  function(x,y,n,j0,e10,e20){
    
    th <- matrix(0,100,5)
    
    # Iteration 0
    
    th[1,1] <- e10
    th[1,2] <- e20
    bc <- beta.calc.E2E2.D(x,y,n,j0,e10,e20)
    th[1,3:5] <- bc$B
    
    # Iterate to convergence (100 Iter max)
    
    for (iter in 2:100){
      m <- iter-1
      
      ec <- eta.calc.D.E2E2(x,y,n,j0,th[m,3:5])
      th[iter,1] <- ec$eta1
      th[iter,2] <- ec$eta2
      
      bc <- beta.calc.E2E2.D(x,y,n,j0,ec$eta1,ec$eta2)
      
      th[iter,3:5] <- bc$B
      theta <- th[1:iter,]
      #delta <- abs(th[iter,]-th[m,])
      delta <- abs(th[iter,]-th[m,])/th[m,]
      if( (delta[1]<.001) & (delta[2]<.001) & (delta[3]<.001)
          & (delta[4]<.001) & (delta[5]<.001) )
        break
    }
    list(theta=theta)
  }

con.est.D.E2E2 <-
  function(xj, est)
  {
    eta1 <- est[1]
    eta2 <- est[2]
    a0 <- est[3]
    a1 <- est[4]
    b1 <- est[5]
    b0 <- a0 + (a1 - b1) * log(xj)
    list(eta1 = eta1, eta2 = eta2, a0 = a0, a1 = a1, b0 = b0, b1 = b1)
  }

con.vals.E2E2.D <-
  function(x, y, n, j, k)
  {
    a <- con.parms.D.E2E2(x, y, n, j, 1, 1)
    nr <- nrow(a$theta)
    est <- a$theta[nr,  ]
    b <- con.est.D.E2E2(x[j], est)
    eta <- c(b$eta1, b$eta2)
    beta <- c(b$a0, b$a1, b$b0, b$b1)
    tau <- x[j]
    list(eta = eta, beta = beta, tau = tau)
  }

p.ll.D <-function(n, j, s2, t2){
  q1 <- n * log(sqrt(2 * pi))
  q2 <- 0.5 * j * (1 + log(s2))
  q3 <- 0.5 * (n - j) * (1 + log(t2))
  - (q1 + q2 + q3)
}

findmax <-function(a){
  maxa<-max(a)
  imax<- which(a==max(a),arr.ind=TRUE)[1]
  jmax<-which(a==max(a),arr.ind=TRUE)[2]
  list(imax = imax, jmax = jmax, value = maxa)
}

beta.calc.E2E2.D <-
  function(x, y, n, j, e1, e2)
  {
    aa <- wmat.E2E2.D(x, y, n, j, e1, e2)
    W <- aa$w
    bb <- rvec.E2E2.D(x, y, n, j, e1, e2)
    R <- bb$r
    beta <- solve(W, R)
    list(B = beta)
  }

eta.calc.D.E2E2 <-
  function(x, y, n, j, theta)
  {
    jp1 <- j + 1
    a0 <- theta[1]
    a1 <- theta[2]
    b1 <- theta[3]
    b0 <- a0 + (a1 - b1) * log(x[j])
    rss1 <- sum((y[1:j] - a0 - a1 * log(x[1:j]))^2)
    rss2 <- sum((y[jp1:n] - b0 - b1 * log(x[jp1:n]))^2)
    e1 <- j/rss1
    e2 <- (n - j)/rss2
    list(eta1 = e1, eta2 = e2)
  }

wmat.E2E2.D <- 
  function(x, y, n, j, e1, e2)
  {
    W <- matrix(0, 3, 3)
    jp1 <- j + 1
    W[1, 1] <- e1 * j + e2 * (n - j) 
    W[1, 2] <- e1 * sum(log(x[1:j])) + e2 * (n - j) * log(x[j])
    W[1, 3] <- e2 * sum(log(x[jp1:n]) - log(x[j]))
    
    W[2, 2] <- e1 * sum(log(x[1:j])^2) + e2 * (n - j) * log(x[j])^2
    W[2, 3] <- e2 * log(x[j]) * sum(log(x[jp1:n]) - log(x[j]))
    
    W[3, 3] <-e2 * sum((log(x[jp1:n]) - log(x[j]))^2)
    
    W[2, 1] <- W[1, 2]
    W[3, 1] <- W[1, 3]
    W[3, 2] <- W[2, 3]
    
    list(w = W)
  }

rvec.E2E2.D <-
  function(x, y, n, j, e1, e2)
  {
    R <- array(0, 3)
    jp1 <- j + 1
    y1j <- sum(y[1:j])
    yjn <- sum(y[jp1:n])
    xy1j <- sum(log(x[1:j]) * y[1:j])
    xyjn <- sum(log(x[jp1:n]) * y[jp1:n])
    
    R[1] <- e1 * y1j + e2 * yjn
    R[2] <- e1 * xy1j + e2 * log(x[j]) * yjn
    R[3] <- e2 * (xyjn - yjn * log(x[j]))
    list(r = R)
  }