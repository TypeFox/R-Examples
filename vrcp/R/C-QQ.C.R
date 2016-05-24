# Nonlinearizable C-QQ
# Common variance for both segments

llsearch.QQ.CC <- function(x, y, n, jlo, jhi)
{
  fj <- matrix(0, n)
  fxy <- matrix(0, jhi - jlo + 1)
  
  jgrid <- expand.grid(jlo:jhi)
  k.ll <- apply(jgrid, 1, p.estFUN.QQ.CC, x = x, y = y, n = n)
  
  fxy <- matrix(k.ll, nrow = jhi-jlo+1)
  rownames(fxy) <- jlo:jhi
  
  z <- findmax(fxy)
  jcrit <- z$imax + jlo - 1
  list(jhat = jcrit, value = max(fxy))
}

#  Function for deriving the ML estimates of the change-points problem.

p.estFUN.QQ.CC <- function(j, x, y, n){
  a <- p.est.QQ.CC(x,y,n,j)
  s2 <- a$sigma2
  return(p.ll.CC(n, j, s2))
}

p.est.QQ.CC <- function(x,y,n,j){
  xa <- x[1:j]
  ya <- y[1:j]
  jp1 <- j+1
  xb <- x[jp1:n]
  yb <- y[jp1:n]
  g1 <- lm(ya ~ xa + I(xa^2))
  ypred <- g1$coef[1]+x[j]*g1$coef[2]+x[j]^2*g1$coef[3]
  ybc <- yb - ypred
  cp <- x[j]
  g2 <- lm(ybc ~ I(xb-cp) + I(xb^2-cp^2))
  b0 <- ypred - g2$coef[2]*cp - g2$coef[3]*cp^2
  beta <-c(g1$coef[1],g1$coef[2],g1$coef[3],b0,g2$coef[2],g2$coef[3])
  s2<- (sum((ya-g1$fit)^2)+sum((yb-g2$fit)^2))/n
  list(a0=beta[1],a1=beta[2],a2=beta[3],b0=beta[4],b1=beta[5],b2=beta[6],sigma2=s2,xj=x[j])
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

