###################################################
### chunk number 1: 
###################################################
#line 6 "cap3.Rnw"
options(prompt="R> ")
options(width=80)


###################################################
### chunk number 2: 
###################################################
#line 31 "cap3.Rnw"
# fig processo 1
set.seed(123)
y <- matrix(, 3, 100)
for(i in 1:3)
 y[i,] <- cumsum(rnorm(100))

y[3,] <- 10+y[3,] 
y[2,] <- -2+y[2,] 
par(mar=c(2,0,0,0))
plot(1:100,y[1,],type="l",axes=FALSE, ylim=c(min(y), max(y)),ylab="",xlab="",xlim=c(0,110))

lines(1:100, y[2,])
lines(1:100, y[3,])

axis(1, c(0,50,150), c(0,expression(bar(t)),""))
text(102, y[1,100], expression(X(t,bar(omega)[1])),adj = c(0, 0.5))
text(102, y[2,100], expression(X(t,bar(omega)[2])),adj = c(0, 0.5))
text(102, y[3,100], expression(X(t,bar(omega)[3])),adj = c(0, 0.5))
abline(v=50, lty=3)
points(c(50,50,50), y[1:3,50])


###################################################
### chunk number 3: 
###################################################
#line 113 "cap3.Rnw"
#pdf("filt1.pdf",width=4,height=2,pointsize=8) 
par(mar=c(0,0,0,0))
plot(c(0,100),c(20,60),type="n", axes=FALSE)#,axes=FALSE,ylab="",xlab="")
lines(c(0,100), c(30,30))
lines(c(50,50), c(29,31))
lines(c(0,0), c(29,31))
lines(c(100,100), c(29,31))
text(20,25, expression(X[1]==0))
text(80,25, expression(X[1]==1))
text(0,35, expression(0))
text(50,35, expression(over(1,2)))
text(100,35, expression(1))

lines(c(0,100), c(50,50))
lines(c(0,0), c(49,51))
lines(c(50,50), c(49,51))
lines(c(25,25), c(49,51))
lines(c(75,75), c(49,51))
lines(c(100,100), c(49,51))
text(10,45, expression(X[2]==0))
text(40,45, expression(X[2]==1))
text(60,45, expression(X[2]==0))
text(90,45, expression(X[2]==1))
text(0,55, expression(0))
text(25,55, expression(over(1,4)))
text(50,55, expression(over(1,2)))
text(75,55, expression(over(3,4)))
text(100,55, expression(1))
#dev.off()


###################################################
### chunk number 4: myfig
###################################################
#line 1135 "cap3.Rnw"
set.seed(123)
n <- 500
t <- 0.3
sim <- 10000
B <- numeric(sim)
for(i in 1:sim){
 X <- sample(c(-1,1), n, replace=TRUE)
 S <- cumsum(X)
 B[i] <- S[n*t]/sqrt(n)
}
plot(density(B),main="",ylab="distribution",
 xlab=expression(B(0.3)),lty=3, axes=FALSE)
axis(1)
g <- function(x) dnorm(x,sd=sqrt(t))
curve( g, -3, 3, add=TRUE)


###################################################
### chunk number 5: 
###################################################
#line 1154 "cap3.Rnw"
par(mar=c(4,3,0,0))
#line 1135 "cap3.Rnw#from line#1155#"
set.seed(123)
n <- 500
t <- 0.3
sim <- 10000
B <- numeric(sim)
for(i in 1:sim){
 X <- sample(c(-1,1), n, replace=TRUE)
 S <- cumsum(X)
 B[i] <- S[n*t]/sqrt(n)
}
plot(density(B),main="",ylab="distribution",
 xlab=expression(B(0.3)),lty=3, axes=FALSE)
axis(1)
g <- function(x) dnorm(x,sd=sqrt(t))
curve( g, -3, 3, add=TRUE)
#line 1156 "cap3.Rnw"


###################################################
### chunk number 6: 
###################################################
#line 1211 "cap3.Rnw"
par(mar=c(4,4,0,0))
rm(list=ls())
tau <- c(20, 30, 10, 90)
n <- length(tau)
N <- 0:n
T <- c(0,cumsum(tau))
plot(T[-(n+1)],N[-(n+1)],axes=FALSE,xlab="t",ylab=expression(N[t]),xlim=c(0,90),ylim=c(0,3),pch=16)

for(i in 2:(n+1)){
 lines( c(T[i-1], T[i]), c(N[i-1],N[i-1]) )
 lines( c(T[i-1], T[i-1]), c(0,N[i-1]), lty=3)
}
points(T[2:n], N[1:(n-1)], pch=1)
axis(1, T, c(0,expression(tau[1]), expression(tau[2]), expression(tau[3]),""))
axis(2, N)


###################################################
### chunk number 7: 
###################################################
#line 1235 "cap3.Rnw"
set.seed(123)
par(mar=c(4,4,0,0))
rm(list=ls())
tau <- c(20, 30, 10, 90)
n <- length(tau)
N <- cumsum(c(0, rnorm(n,sd=2)))
T <- c(0,cumsum(tau))
plot(T[-(n+1)],N[-(n+1)],axes=FALSE,xlab="t",ylab=expression(N[t]),xlim=c(0,90),pch=16)
for(i in 2:(n+1))
 lines( c(T[i-1], T[i]), c(N[i-1],N[i-1]) )
points(T[2:n], N[1:(n-1)], pch=1)
axis(1, T, c(0,expression(tau[1]), expression(tau[2]), expression(tau[3]),""))
axis(2, round(N,2))


