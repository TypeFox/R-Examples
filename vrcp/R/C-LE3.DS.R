# Nonlinearizable C-LE3
# Different variances for two segments
# Smoothness

llsearch.LE3.CDS <- function(x, y, n, jlo, jhi,start_1,start_2,start_3)
{
  fj <- matrix(0, n)
  fxy <- matrix(0, jhi - jlo + 1)
  
  jgrid <- expand.grid(jlo:jhi)
  k.ll <- apply(jgrid, 1, p.estFUN.CDS.LE3, x = x, y = y, n = n,start_1=start_1,start_2=start_2,start_3=start_3)
  
  fxy <- matrix(k.ll, nrow = jhi-jlo+1)
  rownames(fxy) <- jlo:jhi
  
  z <- findmax(fxy)
  jcrit <- z$imax + jlo - 1
  list(jhat = jcrit, value = max(fxy))
}

#  Function for deriving the ML estimates of the change-points problem.

p.estFUN.CDS.LE3 <- function(j, x, y, n,start_1,start_2,start_3){
  a <- p.est.LE3.CDS(x,y,n,j,start_1,start_2,start_3)
  s2 <- a$sigma2
  t2 <- a$tau2
  return(p.ll.CDS(n, j, s2, t2))
}

p.est.LE3.CDS <- function(x,y,n,j,start_1,start_2,start_3){
  xa <- x[1:j]
  ya <- y[1:j]
  jp1 <- j+1
  xb <- x[jp1:n]
  yb <- y[jp1:n]
  fun <- nls(y ~ I(x <= x[j])*(a0 + a1*x) +
               I(x > x[j])*(a0 + a1*x[j] - a1/b2 + a1/b2*exp(b2*(x-x[j]))), 
             start = list(a0 = 10, a1 = -0.8, b2 = 0.5))
  a0 <- summary(fun)$coe[1]
  a1 <- summary(fun)$coe[2]
  b2 <- summary(fun)$coe[3]
  b1 <- a1/b2
  b0 <- a0 + a1 * x[j] - b1
  beta <-c(a0,a1,b0,b1,b2)
  s2<- sum((ya-a0-a1*xa)^2)/j
  t2 <-sum(yb-b0-b1*exp(b2*(xb-x[j]))^2)/(n-j)
  list(a0=beta[1],a1=beta[2],b0=beta[3],b1=beta[4],b2=beta[5],sigma2=s2,tau2=t2,xj=x[j])
}

#  Function to compute the log-likelihood of the change-point problem

p.ll.CDS <- function(n, j, s2, t2){
  q1 <- n * log(sqrt(2 * pi))
  q2 <- 0.5 * j  * (1 + log(s2))
  q3 <- 0.5 * (n - j) * (1 + log(t2))
  - (q1 + q2 + q3)
}

findmax <-function(a)
{
  maxa<-max(a)
  imax<- which(a==max(a),arr.ind=TRUE)[1]
  jmax<-which(a==max(a),arr.ind=TRUE)[2]
  list(imax = imax, jmax = jmax, value = maxa)
}
