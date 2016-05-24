# Method MMP
# Different variances for all three segments 


# mm, nn is the third pair of Gamma parameters
# eta3 ~ Ga(mm,nn)
# Use u and u0 to represent eta3

mmpiter.D <-
function(x,y,n,jlo,jhi,klo,khi,sig0,tau0,u0,rr,ss,vv,ww,mm,nn)
{
xi <- matrix(0,10,5)
xi[1,3] <- sig0
xi[1,4] <- tau0
xi[1,5] <- u0
for (iter in 2:10){
 a <- mmp.D(x,y,n,jlo,jhi,klo,khi,xi[iter-1,3],xi[iter-1,4],xi[iter-1,5],rr,ss,vv,ww,mm,nn)
 xi[iter,1] <- a$jnew
 xi[iter,2] <- a$knew
 xi[iter,3] <- a$signew
 xi[iter,4] <- a$taunew
 xi[iter,5] <- a$unew
 diff1 <- xi[iter,1]-xi[iter-1,1]
 diff2 <- xi[iter,2]-xi[iter-1,2]
 if (abs(diff1) < .1 && abs(diff2) < .1) break
 }
 list(jhat=a$jnew,khat=a$knew,sigma2hat=a$signew,tau2hat=a$taunew,u2hat=a$unew)
}

mmp.D <-
function(x, y, n, jlo, jhi, klo, khi, sig, tau, u, rr, ss, vv, ww, mm, nn)
{
	# 1.  Given sigma^2 & tau^2, find j & k
	step1 <- newjk.D(x, y, n, jlo, jhi, klo, khi, sig, tau, u)
	jhat <- step1$jhat
	khat <- step1$khat

	# 2.  Given j & k, find sigma^2 & tau^2
	step2 <- newst.D(x, y, n, jhat, khat, rr, ss, vv, ww, mm, nn)
	
        sighat <- step2$sig2hat
	tauhat <- step2$tau2hat
        uhat <- step2$u2hat

	list(jnew = jhat, knew = khat, signew = sighat, taunew = tauhat, unew=uhat)
}

newjk.D <-
function(x, y, n, jlo, jhi, klo, khi, sigma2, tau2, u2)
{
	c <- (-10^120)
	fjk <- matrix(c, n, n)
	
	for(j in jlo:jhi) {
		for(k in klo:khi) {
			fjk[j, k] <- lpjk.D(x, y, n, j, k, sigma2, tau2, u2)
		}
	}
	z <- findmax(fjk)
	jcrit <- z$imax
	kcrit <- z$jmax
	list(jhat = jcrit, khat = kcrit, value = max(fjk))
}

newst.D <-
function(x, y, n, j, k, rr, ss, vv, ww, mm, nn)
{
	g <- rss.D(x, y, n, j, k)
	rss1 <- g$S1
	rss2 <- g$S2
	rss3 <- g$S3
	s2 <- (2 * ss + rss1)/(j + 2 * rr + 2)
	t2 <- (2 * ww + rss2)/(k - j + 2 * vv + 2)
  u2 <- (2 * nn + rss3)/(n - k + 2 * mm + 2)
	list(sig2hat = s2, tau2hat = t2, u2hat = u2)
}

lpjk.D <-
function(x, y, n, j, k, sigma2, tau2, u2)
{
	g <- rss.D(x, y, n, j, k)
	rss1 <- g$S1
	rss2 <- g$S2
	rss3 <- g$S3
  ## Yulei's Edit
	x1 <- ((k - j)/2) * log(sigma2/tau2) ## j/2*log(sigma2)+(k-j)/2*log(tau2)-k/2*log(u2)
	x2 <- rss1/(2 * sigma2)
	x3 <- rss2/(2 * tau2)
  x4 <- rss3/(2 * u2)
	x1 - x2 - x3 - x4
}

rss.D <-
function(x, y, n, j, k)
{
	xa <- x[1:j]
	ya <- y[1:j]
	jp1 <- j + 1
	xb <- x[jp1:k]
	yb <- y[jp1:k]
	kp1 <- k + 1
	xc <- x[kp1:n]
	yc <- y[kp1:n]
	g1 <- lm(ya ~ xa)
	g2 <- lm(yb ~ xb)
	g3 <- lm(yc ~ xc)
	beta <- c(g1$coef[1],g1$coef[2],g2$coef[1],g2$coef[2],g3$coef[1],g3$coef[2])
	s1 <- sum((g1$res)^2)
	s2 <- sum((g2$res)^2)
	s3 <- sum((g3$res)^2)
	list(a0=beta[1],a1=beta[2],b0=beta[3],b1=beta[4],c0=beta[5],c1=beta[6], S1 = s1, S2 = s2, 
		S3 = s3)
}

## Yulei's edit to avoid using for-loop
findmax <-function(a){
	maxa<-max(a)
  imax<- which(a==max(a),arr.ind=TRUE)[1]
  jmax<-which(a==max(a),arr.ind=TRUE)[2]
	list(imax = imax, jmax = jmax, value = maxa)
}

