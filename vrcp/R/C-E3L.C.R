# Nonlinearizable C-E3L
# Common variance for both segments

llsearch.E3L.CC <- function(x, y, n, jlo, jhi,start1,start2)
{
  fj <- matrix(0, n)
  fxy <- matrix(0, jhi - jlo + 1)
  
  jgrid <- expand.grid(jlo:jhi)
  k.ll <- apply(jgrid, 1, p.estFUN.E3L.CC, x = x, y = y, n = n,
                start1=start1,start2=start2)
  
  fxy <- matrix(k.ll, nrow = jhi-jlo+1)
  rownames(fxy) <- jlo:jhi
  
  z <- findmax(fxy)
  jcrit <- z$imax + jlo - 1
  list(jhat = jcrit, value = max(fxy))
}

#  Function for deriving the ML estimates of the change-points problem.

p.estFUN.E3L.CC <- function(j, x, y, n, start1,start2){
  a <- p.est.E3L.CC(x,y,n,j,start1,start2)
  s2 <- a$sigma2
  return(p.ll.CC(n, j, s2))
}

p.est.E3L.CC <- function(x,y,n,j,start1,start2){
  xa <- x[1:j]
  ya <- y[1:j]
  jp1 <- j+1
  xb <- x[jp1:n]
  yb <- y[jp1:n] 
  g2 <- lm(yb ~ xb)
  ypred <- g2$coef[1]+x[j]*g2$coef[2]
  changepoint <- x[j]
  fun1<- function(x,a1,a2){ypred - a1 + a1*exp(a2*(x-changepoint))}
  g1 <- nls(ya ~ fun1(xa,a1,a2), data=data.frame(xa,ya),start=list(a1=start1,a2=start2)) 
  beta <-c(summary(g1)$parameter[1],summary(g1)$parameter[2],g2$coef[1],g2$coef[2])
  s2<- (sum((ya-predict(g1, list(x=xa)))^2)+sum((yb-g2$fit)^2))/n
  list(a1=beta[1],a2=beta[2],b0=beta[3],b1=beta[4],sigma2=s2,xj=x[j])
}

#  Function to compute the log-likelihood of the change-point problem

p.ll.CC <- function(n, j, s2){
  q1 <- n * log(sqrt(2 * pi))
  q2 <- 0.5 * n  * (1 + log(s2))
  - (q1 + q2)
}

findmax <-function(a)
{
  maxa<-max(a)
  imax<- which(a==max(a),arr.ind=TRUE)[1]
  jmax<-which(a==max(a),arr.ind=TRUE)[2]
  list(imax = imax, jmax = jmax, value = maxa)
}
