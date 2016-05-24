#  MethUnc-Yulei.s
# Common variance for all three segments


llsearch.C <- function(x, y, n, jlo, jhi, klo, khi,plot)
{
	fjk <- matrix(0, n, n)
	fxy <- matrix(0, (jhi - jlo + 1), (khi - klo + 1))
	
  ## Yulei's edit to avoid using for-loop
	jkgrid <- expand.grid(jlo:jhi, klo:khi)
	k.ll <- apply(jkgrid, 1, p.estFUN.C, x = x, 
	                               y = y, n = n)
	
	fxy <- matrix(k.ll, nrow = jhi-jlo+1, ncol = khi-klo+1)
	rownames(fxy) <- jlo:jhi
	colnames(fxy) <- klo:khi
	if (plot == "TRUE") {
	  jx<-jlo:jhi
	  ky<-klo:khi
	  persp(jx, ky, fxy, xlab = "j", ylab = "k", zlab = "LL(x,y,j,k)")
	  title("Log-likelihood Surface")
	}
	z <- findmax(fxy)
	jcrit <- z$imax + jlo - 1
	kcrit <- z$jmax + klo - 1
	list(jhat = jcrit, khat = kcrit, value = max(fxy))
}

#  S-Plus function for deriving
#  the ML estimates of the 
#  two change-points problem.

p.estFUN.C <- function(jk, x, y, n){
  j = jk[1]
  k = jk[2]
  a <- p.est.C(x,y,n,j,k)
  s2 <- a$sigma2
  return(p.ll.C(n, j, k, s2))
}

p.est.C <- function(x,y,n,j,k){
 xa <- x[1:j]
 ya <- y[1:j]
 jp1 <- j+1
 xb <- x[jp1:k]
 yb <- y[jp1:k]
 kp1 <- k+1
 xc <- x[kp1:n]
 yc <- y[kp1:n]
 g1 <- lm(ya ~ xa)
 g2 <- lm(yb ~ xb)
 g3 <- lm(yc ~ xc)
 beta <-c(g1$coef[1],g1$coef[2],g2$coef[1],g2$coef[2],g3$coef[1],g3$coef[2])
 s2 <- (sum((ya-g1$fit)^2)+sum((yb-g2$fit)^2)+sum((yc-g3$fit)^2))/n
list(a0=beta[1],a1=beta[2],b0=beta[3],b1=beta[4],c0=beta[5],c1=beta[6],sigma2=s2,xj=x[j],xk=x[k])
}

#  S-Plus function to compute the
#  log-likelihood of the two
#  change-point problem

p.ll.C <- function(n, j, k, s2){
 q1 <- n * log(sqrt(2 * pi))
 q2 <- 0.5 * n  * (1 + log(s2))
 - (q1 + q2)
}

##Yulei's edit to avoid using for-loop
findmax <-function(a)
  {
    maxa<-max(a)
    imax<- which(a==max(a),arr.ind=TRUE)[1]
    jmax<-which(a==max(a),arr.ind=TRUE)[2]
	list(imax = imax, jmax = jmax, value = maxa)
}
