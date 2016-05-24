###################################################
### chunk number 1: 
###################################################
#line 6 "cap10.Rnw"
options(prompt="R> ")
options(width=80)


###################################################
### chunk number 2: 
###################################################
#line 144 "cap10.Rnw"
require(sde)
require(fImport)
Delta <- 1/252
S <- yahooSeries("AAPL", from="2009-01-01", to="2009-12-31")
Close <- S[,"AAPL.Close"]
sqrt(var(returns(Close))/Delta)
Close <- rev(Close)
cp <- cpoint(as.ts(Close))
cp


###################################################
### chunk number 3: 
###################################################
#line 163 "cap10.Rnw"
Close1 <- Close[time(Close)[1:cp$k0],]
Close2 <- Close[time(Close)[-(1:cp$k0)],]

sqrt(var(returns(Close1))/Delta)
sqrt(var(returns(Close2))/Delta)


###################################################
### chunk number 4: 
###################################################
#line 171 "cap10.Rnw"
cp2 <- cpoint(as.ts(Close1))
cp2


###################################################
### chunk number 5: cpoint-plot
###################################################
#line 176 "cap10.Rnw"
plot(returns(Close),theme="white",ylab="AAPL Returns",main="",xlab="")
abline(v=time(Close)[cp$k0],lty=3)
abline(v=time(Close1)[cp2$k0],lty=3)


###################################################
### chunk number 6: 
###################################################
#line 183 "cap10.Rnw"
par(mar=c(3,3,1,1))
#line 176 "cap10.Rnw#from line#184#"
plot(returns(Close),theme="white",ylab="AAPL Returns",main="",xlab="")
abline(v=time(Close)[cp$k0],lty=3)
abline(v=time(Close1)[cp2$k0],lty=3)
#line 185 "cap10.Rnw"


###################################################
### chunk number 7: cpoint1
###################################################
#line 404 "cap10.Rnw"
library(yuima)
diff.matrix <- matrix(c("theta1.1*x1","0*x2","0*x1","theta1.2*x2"), 2, 2)
drift.c <- c("1-x1", "3-x2")
drift.matrix <- matrix(drift.c, 2, 1)
ymodel <- setModel(drift=drift.matrix, diffusion=diff.matrix, time.variable="t",
state.variable=c("x1", "x2"), solve.variable=c("x1", "x2"))


###################################################
### chunk number 8: cpoint2 eval=FALSE
###################################################
## #line 412 "cap10.Rnw"
## require(yuima)
## diff.matrix <- matrix(c("theta1.1*x1","0*x2","0*x1","theta1.2*x2"), 2, 2)
## drift.c <- c("1-x1", "3-x2")
## drift.matrix <- matrix(drift.c, 2, 1)
## ymodel <- setModel(drift=drift.matrix, diffusion=diff.matrix, time.variable="t",
## state.variable=c("x1", "x2"), solve.variable=c("x1", "x2"))


###################################################
### chunk number 9: cpoint3
###################################################
#line 422 "cap10.Rnw"
n <- 1000

set.seed(123)

t1 <- list(theta1.1=.1, theta1.2=0.2)
t2 <- list(theta1.1=.6, theta1.2=.6)

tau <- 0.4
ysamp1 <- setSampling(n=tau*n, Initial=0, delta=0.01)
yuima1 <- setYuima(model=ymodel, sampling=ysamp1)
yuima1 <- simulate(yuima1, xinit=c(1, 1), true.parameter=t1)

x1 <- yuima1@data@zoo.data[[1]]
x1 <- as.numeric(x1[length(x1)])
x2 <- yuima1@data@zoo.data[[2]]
x2 <- as.numeric(x2[length(x2)])

ysamp2 <- setSampling(Initial=n*tau*0.01, n=n*(1-tau), delta=0.01)
yuima2 <- setYuima(model=ymodel, sampling=ysamp2)
yuima2 <- simulate(yuima2, xinit=c(x1, x2), true.parameter=t2)
yuima <- yuima1
yuima@data@zoo.data[[1]] <- c(yuima1@data@zoo.data[[1]], yuima2@data@zoo.data[[1]][-1])
yuima@data@zoo.data[[2]] <- c(yuima1@data@zoo.data[[2]], yuima2@data@zoo.data[[2]][-1])


###################################################
### chunk number 10: cpoint4
###################################################
#line 450 "cap10.Rnw"
par(mar=c(3,3,1,1))
plot(yuima)


###################################################
### chunk number 11: cpoint5
###################################################
#line 459 "cap10.Rnw"
t.est <- CPoint(yuima,param1=t1,param2=t2, plot=TRUE)
t.est$tau


###################################################
### chunk number 12: cpoint6
###################################################
#line 464 "cap10.Rnw"
low <- list(theta1.1=0, theta1.2=0)
tmp1 <- qmleL(yuima,start=list(theta1.1=0.3,theta1.2=0.5),t=1.5,lower=low, method="L-BFGS-B")
tmp2 <- qmleR(yuima,start=list(theta1.1=0.3,theta1.2=0.5), t=8.5,lower=low, method="L-BFGS-B")


###################################################
### chunk number 13:  eval=FALSE
###################################################
## #line 469 "cap10.Rnw"
## #line 464 "cap10.Rnw#from line#469#"
## low <- list(theta1.1=0, theta1.2=0)
## tmp1 <- qmleL(yuima,start=list(theta1.1=0.3,theta1.2=0.5),t=1.5,lower=low, method="L-BFGS-B")
## tmp2 <- qmleR(yuima,start=list(theta1.1=0.3,theta1.2=0.5), t=8.5,lower=low, method="L-BFGS-B")
## #line 470 "cap10.Rnw"


###################################################
### chunk number 14: 
###################################################
#line 472 "cap10.Rnw"
coef(tmp1)
coef(tmp2)


###################################################
### chunk number 15: 
###################################################
#line 477 "cap10.Rnw"
t.est1 <- CPoint(yuima,param1=coef(tmp1),param2=coef(tmp2))
t.est1$tau


###################################################
### chunk number 16: 
###################################################
#line 482 "cap10.Rnw"
tmp11 <- qmleL(yuima,start=as.list(coef(tmp1)), t=t.est1$tau-0.1,lower=low, method="L-BFGS-B")
coef(tmp11)
tmp21 <- qmleR(yuima,start=as.list(coef(tmp2)), t=t.est1$tau+0.1,lower=low, method="L-BFGS-B")
coef(tmp21)


###################################################
### chunk number 17: 
###################################################
#line 489 "cap10.Rnw"
t.est2 <- CPoint(yuima,param1=coef(tmp11),param2=coef(tmp21))
t.est2$tau


###################################################
### chunk number 18: 
###################################################
#line 529 "cap10.Rnw"
par(mar=c(0,0,0,0))
plot(1, xlim=c(-1,10), ylim=c(-1,6), type="n",axes=F,ylab="",xlab="")

rect(0, 0, 10, 5, lty=3)
rect(0, 0, 10, 2.5, lty=3)
lines(c(0.5,0.5), c(0,2.5), lty=2)
lines(c(3,3), c(0,2.5), lty=2)
lines(c(5,5), c(0,2.5), lty=2)
lines(c(9,9), c(0,2.5), lty=2)

lines(c(2,2), c(2.5,5), lty=2)
lines(c(4,4), c(2.5,5), lty=2)
lines(c(7,7), c(2.5,5), lty=2)
lines(c(9.5,9.5), c(2.5,5), lty=2)

rect( 4, 2.5, 7, 5, density=20, border=NA)
rect( 3, 0, 5, 2.5, density=20, angle=135, border=NA)
rect( 5, 0, 9, 2.5, density=20, angle=135, border=NA)

text(-1,4.5, expression(X^1))
text(-1,2, expression(X^2))

text(0,-0.3, expression(0))
text(10,-0.3, expression(T))

text(1, 6, expression(I^1))
text(3, 6, expression(I^2))
text(5.5, 6, expression(I^3))
text(8, 6, expression(I^4))
text(9.9, 6, expression(I^5))

text(0.2, -1, expression(J^1))
text(2, -1, expression(J^2))
text(4, -1, expression(J^3))
text(7, -1, expression(J^4))
text(9.5, -1, expression(J^5))

text(2.1, 5.2, expression(T^{11}))
text(4.1, 5.2, expression(T^{12}))
text(7.1, 5.2, expression(T^{13}))
text(9.6, 5.2, expression(T^{14}))

text(0.6, -0.3, expression(T^{21}))
text(3.1, -0.3, expression(T^{22}))
text(5.1, -0.3, expression(T^{23}))
text(9.1, -0.3, expression(T^{24}))


###################################################
### chunk number 19: 
###################################################
#line 646 "cap10.Rnw"
# diffusion coefficient for process 1
diff1 <- function(t,x1=0, x2=0) sqrt(1+t)
# diffusion coefficient for process 2
diff2 <- function(t,x1=0, x2=0) sqrt(1+t^2)
# correlation
rho <- function(t,x1=0, x2=0) sqrt(1/2)
# coefficient matrix for diffusion term
diff.matrix <- matrix( c( "diff1(t,x1,x2)", 
"diff2(t,x1,x2) * rho(t,x1,x2)", "", 
"diff2(t,x1,x2) * sqrt(1-rho(t,x1,x2)^2)"),2,2)
# Model sde using yuima.model
cor.mod <- setModel(drift = c("",""), 
diffusion = diff.matrix, 
solve.variable=c("x1","x2"))


###################################################
### chunk number 20: 
###################################################
#line 663 "cap10.Rnw"
true.theta <- function( T, sigma1, sigma2, rho)
{ 	f <- function(t)  { sigma1(t) * sigma2(t) * rho(t) }
 	integrate(f,0,T)
}


###################################################
### chunk number 21: 
###################################################
#line 695 "cap10.Rnw"
set.seed(123)
Terminal <- 1
n <- 1000
yuima.samp <- setSampling(Terminal=Terminal,n=n)
yuima <- setYuima(model=cor.mod, sampling=yuima.samp)
X <- simulate(yuima)
# Cumulative Covariance
theta <- true.theta(T=Terminal, sigma1=diff1, sigma2=diff2, rho=rho)$value
theta


###################################################
### chunk number 22: 
###################################################
#line 707 "cap10.Rnw"
cce(X)$covmat[1,2]


###################################################
### chunk number 23: 
###################################################
#line 712 "cap10.Rnw"
p1 <- 0.2
p2 <- 0.3
newsamp <- setSampling(
 random=list(rdist=c( function(x) rexp(x, rate=p1*n/Terminal), 
  function(x) rexp(x, rate=p1*n/Terminal))) )


###################################################
### chunk number 24: 
###################################################
#line 720 "cap10.Rnw"
Y <- subsampling(X, sampling=newsamp)


###################################################
### chunk number 25: 
###################################################
#line 724 "cap10.Rnw"
cce(Y)$covmat[1,2]


###################################################
### chunk number 26: 
###################################################
#line 730 "cap10.Rnw"
par(mar=c(0,0,0,0))
plot(Y, oma = c(2, 0, 0.1, 0.1),mar=c(2,4,0,0))


###################################################
### chunk number 27: 
###################################################
#line 739 "cap10.Rnw"
var.c <- function(T, p1,p2, sigma1, sigma2, rho)
{
  tmp_integrand1 <- function(t) (sigma1(t) * sigma2(t))^2
  i1 <- integrate(tmp_integrand1,0,T)
  tmp_integrand2 <- function(t) (sigma1(t) * sigma2(t) * rho(t))^2
  i2 <- integrate(tmp_integrand2,0,T)
  2*(1/p1 + 1/p2)* i1$value + 2*(1/p1+1/p2 - 1/(p1+p2)) * i2$value
}


###################################################
### chunk number 28: 
###################################################
#line 750 "cap10.Rnw"
vc <- var.c(T=Terminal, p1, p2, diff1, diff2, rho)
sqrt(vc/n)


###################################################
### chunk number 29: interest eval=FALSE
###################################################
## #line 868 "cap10.Rnw"
## library(Ecdat)
## library(sde)
## data(Irates)
## rates <- Irates[,"r1"]
## plot(rates)


###################################################
### chunk number 30: 
###################################################
#line 877 "cap10.Rnw"
par(mar=c(3,3,1,1))
#line 868 "cap10.Rnw#from line#878#"
library(Ecdat)
library(sde)
data(Irates)
rates <- Irates[,"r1"]
plot(rates)
#line 879 "cap10.Rnw"


###################################################
### chunk number 31: 
###################################################
#line 890 "cap10.Rnw"
require(yuima)
X <- window(rates, start=1964.471, end=1989.333)
mod <- setModel(drift="alpha+beta*x", diffusion=matrix("sigma*x^gamma",1,1))
yuima <- setYuima(data=setData(X), model=mod)


###################################################
### chunk number 32:  eval=FALSE
###################################################
## #line 896 "cap10.Rnw"
## require(yuima)
## X <- window(rates, start=1964.471, end=1989.333)
## 
## mod <- setModel(drift="alpha+beta*x", diffusion=matrix("sigma*x^gamma",1,1))
## yuima <- setYuima(data=setData(X), model=mod)


###################################################
### chunk number 33: 
###################################################
#line 905 "cap10.Rnw"
lambda1 <- list(alpha=1,beta=1, sigma =1, gamma =1)
start <- list(alpha=1, beta =-.1, sigma =.1, gamma =1)
low <- list(alpha=-5, beta =-5, sigma =-5, gamma =-5)
upp <- list(alpha=8, beta =8, sigma =8, gamma =8)
lasso1 <- lasso(yuima, lambda1, start=start, lower=low, upper=upp, method="L-BFGS-B")


###################################################
### chunk number 34: 
###################################################
#line 913 "cap10.Rnw"
round(lasso1$mle, 3)


###################################################
### chunk number 35: 
###################################################
#line 917 "cap10.Rnw"
round(lasso1$lasso, 3)


###################################################
### chunk number 36: 
###################################################
#line 921 "cap10.Rnw"
lambda10 <- list(alpha=10, beta =10, sigma =10, gamma =10)
lasso10 <- lasso(yuima, lambda10, start=start, lower=low, upper=upp, method="L-BFGS-B")


###################################################
### chunk number 37: 
###################################################
#line 926 "cap10.Rnw"
round(lasso10$mle, 3)
round(lasso10$lasso, 3)


###################################################
### chunk number 38: 
###################################################
#line 935 "cap10.Rnw"
theta.tilde <- sprintf("%.4f",lasso1$mle[c("alpha","beta","sigma", "gamma")] )
theta.hat1 <- sprintf("%.4f",lasso1$lasso[c("alpha","beta","sigma", "gamma")])
theta.hat2 <- sprintf("%.4f",lasso10$lasso[c("alpha","beta","sigma", "gamma")])
SIGMA <- sprintf("%.3f",lasso1$sd.mle[c("alpha","beta","sigma", "gamma")] )
SIGMA1 <- sprintf("%.3f",lasso1$sd.lasso[c("alpha","beta","sigma", "gamma")] )
SIGMA2 <- sprintf("%.3f",lasso10$sd.lasso[c("alpha","beta","sigma", "gamma")] )


###################################################
### chunk number 39: multidim
###################################################
#line 993 "cap10.Rnw"
diff.matrix <- matrix(c("theta1.1","theta1.2", "1", "1"), 2, 2)

drift.c <- c("-theta2.1*x1", "-theta2.2*x2", "-theta2.2", "-theta2.1")
drift.matrix <- matrix(drift.c, 2, 2)

ymodel <- setModel(drift=drift.matrix, diffusion=diff.matrix, time.variable="t",
                   state.variable=c("x1", "x2"), solve.variable=c("x1", "x2"))
n <- 100
ysamp <- setSampling(Terminal=(n)^(1/3), n=n)
yuima <- setYuima(model=ymodel, sampling=ysamp)


###################################################
### chunk number 40: multidimsim
###################################################
#line 1006 "cap10.Rnw"
set.seed(123)
truep <- list(theta1.1=0.6, theta1.2=0,theta2.1=0.5, theta2.2=0)
yuima <- simulate(yuima, xinit=c(1, 1), true.parameter=truep)


###################################################
### chunk number 41: lasso
###################################################
#line 1012 "cap10.Rnw"
est <- lasso(yuima, start=list(theta2.1=0.8, theta2.2=0.2, theta1.1=0.7, theta1.2=0.1),
 lower=list(theta1.1=1e-10,theta1.2=1e-10,theta2.1=.1,theta2.2=1e-10),
 upper=list(theta1.1=4,theta1.2=4,theta2.1=4,theta2.2=4), method="L-BFGS-B")

# TRUE
unlist(truep)

# QMLE
round(est$mle,3)

# LASSO
round(est$lasso,3) 


###################################################
### chunk number 42: 
###################################################
#line 1122 "cap10.Rnw"
require(sde)
data(quotes)
Series <- quotes
nSeries <- dim(Series)[2]
plot(Series,main="",xlab="")


###################################################
### chunk number 43: tseries eval=FALSE
###################################################
## #line 1129 "cap10.Rnw"
## plot(Series,main="",xlab="")


###################################################
### chunk number 44: 
###################################################
#line 1134 "cap10.Rnw"
par(mar=c(0,0,0,0))
#line 1129 "cap10.Rnw#from line#1135#"
plot(Series,main="",xlab="")
#line 1136 "cap10.Rnw"


###################################################
### chunk number 45: 
###################################################
#line 1143 "cap10.Rnw"
STSdist <- function(data){
 nSer <-NCOL(data)
 d <- matrix(0, nSer, nSer) 
 colnames(d) <- colnames(data)
 rownames(d) <- colnames(data)
 DELTA <- deltat(data)
 for(i in 1:(nSer-1))
  for(j in (i+1):nSer){
   d[i,j] <- sqrt(sum((diff(data[,i])/DELTA-diff(data[,j])/DELTA)^2))  
   d[j,i] <- d[i,j]
  }
 invisible(d)
}


###################################################
### chunk number 46: MOdist
###################################################
#line 1160 "cap10.Rnw"
dMO <- MOdist(Series)
dMO <- dMO/max(dMO)

dEUC <- dist(t(Series))
dEUC <- dEUC/max(dEUC)

dSTS <- STSdist(Series)
dSTS <- dSTS/max(dSTS)

require(dtw)
dDTW <- dist(t(Series), method="dtw")
dDTW <- dDTW/max(dDTW)


###################################################
### chunk number 47: dendro2 eval=FALSE
###################################################
## #line 1175 "cap10.Rnw"
## cl <- hclust(dMO) 
## plot(cl,main="Markov Operator Distance",xlab="",ylim=c(0,1))
## rect.hclust(cl,k=6,border=gray(0.5))
## 
## cl1 <- hclust(as.dist(dEUC)) 
## plot(cl1, main="Euclidean Distance",xlab="",ylim=c(0,1))
## rect.hclust(cl1,k=6,border=gray(0.5))
## 
## 
## cl2 <- hclust(as.dist(dSTS)) 
## plot(cl2, main="STS Distance",xlab="",ylim=c(0,1))
## rect.hclust(cl2,k=6,border=gray(0.5))
## 
## cl3 <- hclust(as.dist(dDTW)) 
## plot(cl3, main="DTW Distance",xlab="",ylim=c(0,1))
## rect.hclust(cl3,k=6,border=gray(0.5))


###################################################
### chunk number 48: 
###################################################
#line 1200 "cap10.Rnw"
par(mfrow=c(2,2))
par(mar=c(1,3,3,0))
#line 1175 "cap10.Rnw#from line#1202#"
cl <- hclust(dMO) 
plot(cl,main="Markov Operator Distance",xlab="",ylim=c(0,1))
rect.hclust(cl,k=6,border=gray(0.5))

cl1 <- hclust(as.dist(dEUC)) 
plot(cl1, main="Euclidean Distance",xlab="",ylim=c(0,1))
rect.hclust(cl1,k=6,border=gray(0.5))


cl2 <- hclust(as.dist(dSTS)) 
plot(cl2, main="STS Distance",xlab="",ylim=c(0,1))
rect.hclust(cl2,k=6,border=gray(0.5))

cl3 <- hclust(as.dist(dDTW)) 
plot(cl3, main="DTW Distance",xlab="",ylim=c(0,1))
rect.hclust(cl3,k=6,border=gray(0.5))
#line 1203 "cap10.Rnw"


###################################################
### chunk number 49: 
###################################################
#line 1223 "cap10.Rnw"
Sim <- function(g1, g2){
	G1 <- unique(g1)
	G2 <- unique(g2)
	l1 <- length(G1)
	l2 <- length(G2)	
	sim <- matrix(, l1, l2)
	for(i in 1:l1){
		idx <- which(g1==i)
		for(j in 1:l2){
			idx2 <- which(g2==j)
			sim[i,j] <- 2*length(intersect(idx, idx2))/(length(idx)+length(idx2))
		}
	}
	sum(apply(sim, 2, max))/l1
}
G <- cutree(cl,k=6)
G1 <- cutree(cl1,k=6)
G2 <- cutree(cl2,k=6)
G3 <- cutree(cl3,k=6)


A <- matrix(,4,4)
A[1,1] <- Sim(G,G)
A[1,2] <- Sim(G,G1)
A[1,3] <- Sim(G,G2)
A[1,4] <- Sim(G,G3)
A[2,1] <- Sim(G1,G)
A[2,2] <- Sim(G1,G1)
A[2,3] <- Sim(G1,G2)
A[2,4] <- Sim(G1,G3)
A[3,1] <- Sim(G2,G)
A[3,2] <- Sim(G2,G1)
A[3,3] <- Sim(G2,G2)
A[3,4] <- Sim(G2,G3)
A[4,1] <- Sim(G3,G)
A[4,2] <- Sim(G3,G1)
A[4,3] <- Sim(G3,G2)
A[4,4] <- Sim(G3,G3)

S <- (A+t(A))/2


###################################################
### chunk number 50: 
###################################################
#line 1266 "cap10.Rnw"
A <- round(A,2)
S <-round(S,2)


###################################################
### chunk number 51: tseries2 eval=FALSE
###################################################
## #line 1289 "cap10.Rnw"
## plot(Series,ylim=range(Series), main="",xlab="")


###################################################
### chunk number 52: 
###################################################
#line 1294 "cap10.Rnw"
par(mar=c(0,0,0,0))
#line 1289 "cap10.Rnw#from line#1295#"
plot(Series,ylim=range(Series), main="",xlab="")
#line 1296 "cap10.Rnw"


