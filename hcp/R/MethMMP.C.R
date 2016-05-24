#  MethMMP.s  
#
#  S-Plus functions used in the
#  Maximization-Maximization-Posterior method

#  mmpiter.s
mmpiter.C <- function(x,y,n,jlo,jhi,klo,khi,sig0,rr,ss){
xi <- matrix(0,10,4)
xi[1,3] <- sig0

for (iter in 2:10){
  a <- mmp.C(x,y,n,jlo,jhi,klo,khi,xi[iter-1,3],rr,ss)
  xi[iter,1] <- a$jnew
  xi[iter,2] <- a$knew
  xi[iter,3] <- a$signew
  diff1 <- xi[iter,1]-xi[iter-1,1]
  diff2 <- xi[iter,2]-xi[iter-1,2]
  if (abs(diff1) < .1 && abs(diff2) < .1) break
}
list(jhat=a$jnew,khat=a$knew,sigma2hat=a$signew)
}

#  mmp.s
mmp.C <- function(x, y, n, jlo, jhi, klo, khi, sig, rr, ss)
  {
    
	# 1.  Given sigma^2 & tau^2, find j & k
	step1 <- newjk.C(x, y, n, jlo, jhi, klo, khi, sig)

	# 2.  Given j & k, find sigma^2 & tau^2
	step2 <- newst.C(x, y, n, step1$jhat, step1$khat, rr, ss)

	list(jnew = step1$jhat, knew = step1$khat, signew = step2$sig2hat)
}

#  newjk.s:  Maximize log[pi(j,k|...)]
newjk.C <- function(x, y, n, jlo, jhi, klo, khi, sigma2)
{
	fjk <- matrix(-10^120, n, n)
	
	## Yulei's edit to avoid using for-loop
	jkgrid <- expand.grid(jlo:jhi, klo:khi)
  res <- data.frame(k.ll = apply(jkgrid, 1, lpjk.C, x = x, 
	                               y = y, n = n,sigma2=sigma2))
	
	fjk[jlo:jhi,klo:khi] <- matrix(res$k.ll, nrow = jhi-jlo+1, ncol = khi-klo+1)	
	z<-findmax(fjk)
	list(jhat = z$imax, khat = z$jmax, value = max(fjk))
}

#  newst.s:  Maximize pi(sigma^2,tau^2|...)
newst.C <- function(x, y, n, j0, k0, rr, ss)
{
	g <- rss.C(x,y,n,j0,k0)
	s2 <- (2 * ss + g$S1 + g$S2 + g$S3)/(n + 2 * rr + 2)
	list(sig2hat = s2)
}

#  log p_jk
##Yulei's edit
lpjk.C <-
  function(jk, x, y, n, sigma2)
  {  
    j<-jk[1]
    k<-jk[2]
    g <- rss.C(x, y, n, j, k)
    x1 <- (n/2)* log(sigma2)
    x2 <- (g$S1 + g$S2 + g$S3)/(2 * sigma2)
    x1 - x2
  }

#  rss - Residual Sum of Squares
rss.C <- function(x,y,n,j0,k0)
{
	jp1 <- j0 + 1
	kp1 <- k0 + 1 
	g1 <- lm(y[1:j0] ~ x[1:j0])
	g2 <- lm(y[jp1:k0] ~ x[jp1:k0])
	g3 <- lm(y[kp1:n] ~ x[kp1:n])
	s1 <- sum((g1$res)^2)
	s2 <- sum((g2$res)^2)
	s3 <- sum((g3$res)^2)
	list(S1 = s1, S2 = s2, S3 = s3)
}

## Yulei's edit to avoid using for-loop
findmax <-function(a){
	maxa<-max(a)
  imax<- which(a==max(a),arr.ind=TRUE)[1]
  jmax<-which(a==max(a),arr.ind=TRUE)[2]
	list(imax = imax, jmax = jmax, value = maxa)
} 