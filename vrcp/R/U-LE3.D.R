# Linearizable U-LE3
# Different variances for two segments

llsearch.LE3.D <- function(x, y, n, jlo, jhi,start1, start2, start3)
{
  fj <- matrix(0, n)
  fxy <- matrix(0, jhi - jlo + 1)
  
  jgrid <- expand.grid(jlo:jhi)
  k.ll <- apply(jgrid, 1, p.estFUN.LE3.D, x = x, y = y, n = n,
  start1 = start1, start2=start2, start3=start3)
  
  fxy <- matrix(k.ll, nrow = jhi-jlo+1)
  rownames(fxy) <- jlo:jhi
  
  z <- findmax(fxy)
  jcrit <- z$imax + jlo - 1
  list(jhat = jcrit, value = max(fxy))
}

#  Function for deriving the ML estimates of the change-points problem.

p.estFUN.LE3.D <- function(j, x, y, n,start1, start2, start3){
  a <- p.est.LE3.D(x,y,n,j,start1, start2, start3)
  s2 <- a$sigma2
  t2 <- a$tau2
  return(p.ll.D(n, j, s2, t2))
}

p.est.LE3.D <- function(x,y,n,j,start1, start2, start3){
  xa <- x[1:j]
  ya <- y[1:j]
  jp1 <- j+1
  xb <- x[jp1:n]
  yb <- y[jp1:n]
  g1 <- lm(ya ~ xa)
  fun2<-function(x,b0,b1,b2){b0*(1-exp(b1*(x-b2)))}
  g2=nls(yb~fun2(xb,b0,b1,b2),data=data.frame(xb,yb),
         start=list(b0=start1,b1=start2,b2=start3))
  beta <-c(g1$coef[1],g1$coef[2],summary(g2)$parameter[1],summary(g2)$parameter[2],
           summary(g2)$parameter[3])
  s2<- sum((ya-g1$fit)^2)/j
  t2 <-sum((yb-predict(g2,list(x=xb)))^2)/(n-j)
  list(a0=beta[1],a1=beta[2],b0=beta[3],b1=beta[4],b2=beta[5],sigma2=s2,tau2=t2,xj=x[j])
}

#  Function to compute the log-likelihood of the change-point problem

p.ll.D <- function(n, j, s2, t2){
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





