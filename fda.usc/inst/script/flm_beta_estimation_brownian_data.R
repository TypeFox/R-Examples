################################################################################
# Estimation functional beta parameter  using simulated data
# where the theoretical  beta is known
#
# For more details see: Febrero-Bande, M., Galeano, P. and Gonzalez-Manteiga, W.
# (2010). Measures of influence for the functional linear model with scalar
# response. Journal of Multivariate Analysis 101, 327-339.
################################################################################


#####################################################################################
# The function fdata.brown() generates functional data independent of a movement
# Brownian defined in the interval [0, Tend]. The data appear in a matrix
# (J x n) size where n is the number of curves and J is the number of discretized
# points, which is supposed to be the same for all curves.
#####################################################################################
fdata.brown <- function(n=1000,J=101,tt=seq(0,1,len=J),rtt=range(tt),names=NULL){
	X <- matrix(NA,ncol=J,nrow=n)
	for (i in 1:n){
     		X[i,1] <- 0
      	for (j in 2:J){X[i,j] <- X[i,j-1] + rnorm(1,0,sqrt(tt[j]-tt[j-1]))}
	}
	X<-fdata(X,tt,rtt,names=names)
	return(X)
}
#####################################################################################
# The function fdata.cfs.2003.a()   generates simulated data of the functional
# linear linear model with scalar response, Cardot, Ferraty y Sarda (2003).
#####################################################################################
 fdata.cfs.2003.a <- function(x=NULL,R2){
	if  (is.null(x)) x<- fdata.brown()
	tj <- x$argvals
	rtt<-x$rangeval
	X <- x$data
	J<-ncol(x)
	n<-nrow(x)
	v1 <- sqrt(2) * sin(.5 * pi * tj)
	v2 <- sqrt(2) * sin(1.5 * pi * tj)
	v3 <- sqrt(2) * sin(2.5 * pi * tj)
	lamb1 <- 1 / (.5 * pi)^2
	lamb2 <- 1 / (1.5 * pi)^2
	lamb3 <- 1 / (2.5 * pi)^2
	b1 <- 2 / sqrt(2)
	b2 <- 4 / sqrt(2)
	b3 <- 5 / sqrt(2)
	bet <- b1 * v1 + b2 * v2 + b3 * v3
	varXg <- b1^2 * lamb1 + b2^2 * lamb2 + b3^2 * lamb3
	s2eps <- varXg * (1 / R2 - 1)
	Xmean <- fdata.cen(x)
	yp <- (Xmean$Xcen$data %*% bet) * (1 / J)
	e <- matrix(rnorm(n,0,sqrt(s2eps)),ncol=1)
	y <- yp + e
	print("out")
	return(list("x"=x,"y"=y,"v1"=fdata(v1,tj,rtt),"v2"=fdata(v2,tj,rtt),
  "v3"=fdata(v3,tj,rtt),"bet"=fdata(bet,tj,rtt),"s2eps"=s2eps))
}
#####################################################################################
# The function data.cfs.2003.b()  generates simulated data of the functional
# linear linear model with scalar response, Cardot, Ferraty y Sarda (2003).
#####################################################################################
 fdata.cfs.2003.b <- function(x=NULL,R2){
	if  (is.null(x)) x<- fdata.brown()
	tj <- x$argvals
	rtt<-x$rangeval
	X <- x$data
	J<-ncol(x)
	n<-nrow(x)
	maxK <- qr(X)$rank
	varXg <- 0
	bet <- log(15*tj^2 + 10) + cos(4*pi*tj)
	for (k in 1:maxK){
		integrand <- function(t) {sqrt(2)*(log(15*t^2+10)+cos(4*pi*t))*sin((k-.5)*pi*t)}
		bk <- integrate(integrand, lower = rtt[1], upper = rtt[2])$value
		varXg <- varXg + bk^2 / ((k - .5) * pi)^2
	}
	s2eps <- varXg * (1 / R2 - 1)
	Xmean <- fdata.cen(x)
	yp <- (Xmean$Xcen$data %*% bet) * (1 / J)
	e <- matrix(rnorm(n,0,sqrt(s2eps)),ncol=1)
	y <- yp + e
	return(list("x"=x,"y"=y,"bet"=fdata(bet,tj,rtt),"s2eps"=s2eps))
}


#####################################################################################
# The function data.hh.2006 generates simulated data of the functional linear
# model with scalar response, de Hall and Housseini-Nasab (2006).
#####################################################################################
fdata.hh.2006 <- function(n,J,R2){
	tj <- seq(0,1,length.out=J)
	maxK <- min(c(n,J-1))
	v <- matrix(NA,ncol=maxK,nrow=J)
	bet <- pi^2 * (tj^2 - 1/3)
	varXg <- 0
	for (k in 1 : maxK){
		v[,k] <- sqrt(2) * cos(k * pi * tj)
		bk <- (-1)^k * k^(-2) * 2^(3/2)
		varXg <- varXg + bk^2 / k^2
	}
	s2eps <- varXg * (1 / R2 - 1)
  X <- matrix(0,ncol=J,nrow=n)
	for (i in 1 : n){for (k in 1:maxK){X[i,] <- X[i,] + rnorm(1,0,k^(-1)) * t(v[,k])}}
  x<-fdata(X,tj)
  Xmean <- fdata.cen(x)
	yp <- (Xmean$Xcen$data %*% bet) * (1 / J)
	e <- matrix(rnorm(n,0,sqrt(s2eps)),ncol=1)
	y <- yp + e
	R2 <- varXg/(1 + varXg)
	return(list("x"=x,"y"=y,"bet"=fdata(bet,tj),"s2eps"=s2eps))
}

#####################################################################################
# The function data.cfs.2003.c()generates simulated data of the functional linear
# model with scalar response, Cardot, Ferraty y Sarda (2003), but with the
# inclusion of different eigenfunctions.
#####################################################################################
 fdata.cfs.2003.c <- function(x=NULL,R2){
	if  (is.null(x)) x<- fdata.brown()
	tj <- x$argvals
	rtt<-x$rangeval
	X <- x$data
	J<-ncol(x)
	n<-nrow(x)
	v3 <- sqrt(2) * sin(2.5 * pi * tj)
	v5 <- sqrt(2) * sin(4.5 * pi * tj)
	v7 <- sqrt(2) * sin(6.5 * pi * tj)
	lamb3 <- 1 / (2.5 * pi)^2
	lamb5 <- 1 / (4.5 * pi)^2
	lamb7 <- 1 / (6.5 * pi)^2
	b3 <- 2 / sqrt(2)
	b5 <- 4 / sqrt(2)
	b7 <- 5 / sqrt(2)
	bet <- b3 * v3 + b5 * v5 + b7 * v7
	varXg <- b3^2 * lamb3 + b5^2 * lamb5 + b7^2 * lamb7
	s2eps <- varXg * (1 / R2 - 1)
  Xmean <- fdata.cen(X)
	yp <- (Xmean$Xcen$data %*% bet) * (1 / J)
	e <- matrix(rnorm(n,0,sqrt(s2eps)),ncol=1)
	y <- yp + e
	return(list("x"=x,"y"=y,"v3"=v3,"v5"=v5,"v7"=v5,"bet"=fdata(bet,tj,rtt),
  "s2eps"=s2eps))
}
################################################################################
# Case studied
n=100;J=50;r2=0.3
fx=fdata.brown(n,J)
fdatos.a=fdata.cfs.2003.a(fx,R2=r2)
fdatos.b=fdata.cfs.2003.b(fx,R2=r2)
fdatos.c=fdata.cfs.2003.c(fx,R2=r2)
fdatos.d<-fdata.hh.2006(n,J,R2=r2)

################################################################################
# Example for the data generated in Figure 1 Febrero et al. (2010)
par(mfrow=c(2,2))
plot(fdatos.a$bet,main="cfs.2003.a")
plot(fdatos.b$bet,col=2,main="cfs.2003.b")
plot(fdatos.c$bet,col=3,main="cfs.2003.c")
plot(fdatos.d$bet,col=4,main="hh.2006")
par(mfrow=c(2,2))
plot(fdatos.a$x,main="cfs.2003.a")
plot(fdatos.b$x,col=2,main="cfs.2003.b")
plot(fdatos.c$x,col=3,main="cfs.2003.c")
plot(fdatos.d$x,col=4,main="hh.2006")

################################################################################
# Example estimation functional beta parameter  using  simulated data  generates
# via fdata.cfs.2003.a()
obj<-fdatos.a
x<-obj$x
y<-obj$y
dev.off()
plot(obj$bet,main="Beta theoretical",ylab="Beta(t)")

################################################################################
# Functional linear model with principal components
res=fregre.pc(x,y,1:3)
# Displaying Theoretical and estimated  beta
plot(c(obj$bet,res$beta.est),lwd=1:2,main="Beta(t) - Beta.est(t)",ylab="beta(t)")

################################################################################
# n betas estimated via cross-validation
res.infl<-influence.fdata(res) # This take a lot
plot(res.infl$betas,type="l",col=4,main="betas.CV")
lines(res$beta.est,type="l",col=2)
lines(obj$bet,type="l",lwd=2)

################################################################################
# mue.boot betas estimated by smoothing boostrrap

#ARGUMENTS for influence.quan:
#smo:  Smoothing parameter as a proportion of response variance.
#smoX: Smoothing parameter for fdata object as a proportion of variance-covariance
#      matrix of the explanatory functional variable.
#alpha:Significance level.

# This take a lot
resquan=influence.quan(res,res.infl,mue.boot=10, smo = 0.1, smoX = 0.05,
                        alpha=0.95,kmax.fix=TRUE)
plot(resquan$betas.boot,col=4,main="Betas CV, 3 PC, nc=125")
lines(res$beta.est,col=2,lwd=2,lty=2)
lines(obj$bet,type="l",lwd=2)

################################################################################
#  Calculation of the difference in L2-norm  between the true beta
# and estimated by bootstrap
matbeta.est<-obj$bet
for (i in 2:nrow(resquan$betas.boot)) matbeta.est<-c(matbeta.est,obj$bet)
bb<-matbeta.est-resquan$betas.boot
norm.boot<-  norm.fdata(bb)

#  Calculation of the difference in L2-norm  between the beta estimated
#  in fitted model  and the betas estimated by bootstrap
matbeta.est2<-res$beta.est
for (i in 2:nrow(resquan$betas.boot)) matbeta.est2<-c(matbeta.est2,res$beta.est)
bb2<-matbeta.est2-resquan$betas.boot
norm.boot2<-  norm.fdata(bb2)

# displaying the  differencies
plot(density(norm.boot2),main="L2-norm(beta.est-beta.boost)",xlab="DISTANCE")
# percentile 95% of the distances (red line)
dist095<-quantile(norm.boot2,probs=0.95)
abline(v=dist095,col=2)

################################################################################
################################################################################
# The above procedure can be repeated for the function fregre.basis()
# instead of the functon fregre.pc.cv()
################################################################################
################################################################################

################################################################################
# Functional linear model with basis representation
res03=fregre.basis.cv(x,y,11,3,type.basis="fourier")
plot(obj$bet,main="Beta: theoretical and estimated",ylab="Beta(t)")
lines(res03$beta.est,col=3)

################################################################################
# n betas estimated via cross-validation
res03.infl<-influence.fdata(res03)
plot(res03.infl$betas,type="l",col=4,main="Betas CV, 3 Basis of Fourier")
lines(res03$beta.est,type="l",col=3,lty=2)
lines(obj$bet,type="l",lwd=2)

################################################################################
# mue.boot betas estimated by smoothing boostrrap
resquan03=influence.quan(res03,res03.infl,mue.boot=10,kmax.fix=TRUE)
plot(resquan03$betas.boot,col=4,main="Betas bootstrap,3 Bases Fourier")
lines(res03$beta.est,col=2,lwd=2,lty=2)
lines(obj$bet,type="l",lwd=2)


################################################################################
#  Calculation of the difference in L2-norm  between the true beta
# and estimated by bootstrap
plot(density(resquan03$norm.boot),main="norma(beta.est-beta.boost)",xlab="DISTANCE")
dist095<-quantile(resquan03$norm.boot,probs=0.95)
abline(v=dist095,col=2)

################################################################################
# In addition, you can calculate confidence intervals for the basis coefficients
alfa<-0.05
p1=alfa/ncol(resquan03$coefs.boot)
p2=1-alfa/ncol(resquan03$coefs.boot)
IC<-apply(resquan03$coefs.boot,2,quantile,probs=c(p1,p2))
matplot(t(resquan03$coefs.boot),type="l",col=1)
matplot(t(IC),col=2,type="l",add=TRUE)



