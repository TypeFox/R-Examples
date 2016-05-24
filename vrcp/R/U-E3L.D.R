# Linearizable U-E3L
# Different variances for two segments

llsearch.E3L.D <- function(x, y, n, jlo, jhi)
{
  fj <- matrix(0, n)
  fxy <- matrix(0, jhi - jlo + 1)
  
  jgrid <- expand.grid(jlo:jhi)
  k.ll <- apply(jgrid, 1, p.estFUN.E3L.D, x = x, y = y, n = n)
  
  fxy <- matrix(k.ll, nrow = jhi-jlo+1)
  rownames(fxy) <- jlo:jhi
  
  z <- findmax(fxy)
  jcrit <- z$imax + jlo - 1
  list(jhat = jcrit, value = max(fxy))
}

#  Function for deriving the ML estimates of the change-points problem.

p.estFUN.E3L.D <- function(j, x, y, n){
  a <- p.est.E3L.D(x,y,n,j)
  s2 <- a$sigma2
  t2<- a$tau2
  return(p.ll.D(n, j, s2, t2))
}

p.est.E3L.D <- function(x,y,n,j){
  xa <- x[1:j]
  ya <- y[1:j]
  jp1 <- j+1
  xb <- x[jp1:n]
  yb <- y[jp1:n]
  fun1<-function(x,a0,a1,a2){a0*(1-exp(a1*(x-a2)))}
  g1 <- nls(ya~fun1(xa,a0,a1,a2),data=data.frame(xa,ya),start=list(a0=5,a1=0.5,a2=-0.2))
  g2 <- lm(yb ~ xb)
  beta <-c(summary(g1)$parameter[1],summary(g1)$parameter[2],
           summary(g1)$parameter[3],g2$coef[1],g2$coef[2])
  s2<- sum((ya-predict(g1, list(x=xa)))^2)/j
  t2 <- sum((yb-g2$fit)^2)/(n-j)
  list(a0=beta[1],a1=beta[2],a2=beta[3],b0=beta[4],b1=beta[5],sigma2=s2,tau2=t2, xj=x[j])
}

#  Function to compute the log-likelihood of the change-point problem

p.ll.D <- function(n, j, s2, t2){
  q1 <- n * log(sqrt(2 * pi))
  q2 <- 0.5 * n  * (1 + log(s2))
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






