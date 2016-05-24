# Nonlinearizable C-QE3
# Common variances for two segments
# Smoothness

llsearch.QE3.CCS <- function(x, y, n, jlo, jhi,start1,start2,start3,start4)
{
  fj <- matrix(0, n)
  fxy <- matrix(0, jhi - jlo + 1)
  
  jgrid <- expand.grid(jlo:jhi)
  k.ll <- apply(jgrid, 1, p.estFUN.QE3.CCS, x = x, y = y, n = n,start1=start1,start2=start2,start3=start3,start4=start4)
  
  fxy <- matrix(k.ll, nrow = jhi-jlo+1)
  rownames(fxy) <- jlo:jhi
  
  z <- findmax(fxy)
  jcrit <- z$imax + jlo - 1
  list(jhat = jcrit, value = max(fxy))
}

#  Function for deriving the ML estimates of the change-points problem.

p.estFUN.QE3.CCS <- function(j, x, y, n,start1,start2,start3,start4){
  a <- p.est.QE3.CCS(x,y,n,j,start1,start2,start3,start4)
  s2 <- a$sigma2
  return(p.ll.CCS(n, j, s2))
}

p.est.QE3.CCS <- function(x,y,n,j,start1,start2,start3,start4){
  xa <- x[1:j]
  ya <- y[1:j]
  jp1 <- j+1
  xb <- x[jp1:n]
  yb <- y[jp1:n]
  fun <- nls(y ~ I(x <= x[j])*(a0 + a1*x + a2*x^2) +
               I(x > x[j])*(a0 + a1*x[j]+a2*x[j]^2 - (a1+2*a2*x[j])/b2 + (a1+2*a2*x[j])/b2*exp(b2*(x-x[j]))), 
             start = list(a0 = 5, a1 = -2, a2= -1, b2 = 0.3))
  a0 <- summary(fun)$coe[1]
  a1 <- summary(fun)$coe[2]
  a2 <- summary(fun)$coe[3]
  b2 <- summary(fun)$coe[4]  
  b1 <- (a1+2*a2*x[j])/b2
  b0 <- a0 + a1 * x[j] + a2 * x[j]^2 - b1  
  beta <-c(a0, a1, a2, b0, b1, b2)
  s2<- (sum((ya-a0 - a1*xa - a2*xa^2)^2)+sum((yb-b0-b1*exp(b2*(xb-x[j])))^2))/n
  list(a0=beta[1],a1=beta[2],a2=beta[3],b0=beta[4],b1=beta[5],b2=beta[6],sigma2=s2,xj=x[j])
}

#  Function to compute the log-likelihood of the change-point problem

p.ll.CCS <- function(n, j, s2){
  q1 <- n * log(sqrt(2 * pi))
  q2 <- 0.5 * j  * (1 + log(s2))
  - (q1 + q2)
}

findmax <-function(a)
{
  maxa<-max(a)
  imax<- which(a==max(a),arr.ind=TRUE)[1]
  jmax<-which(a==max(a),arr.ind=TRUE)[2]
  list(imax = imax, jmax = jmax, value = maxa)
}