# Chapter 1 R commands

# Section 1.4.2

str(log)
a=c(2,6,-4,9,18)
x <- c(3,6,9)
d=a[c(1,3,5)]
e=3/d
e=lgamma(e^2)

S=readline(prompt="Type  <Return>   to continue : ")

x=matrix(1:4,ncol=3)
print(x[x>5])
print(x[1.])

system.time(crossprod(1:10^6,1:10^6))
system.time(t(1:10^6)%*%(1:10^6))

S=readline(prompt="Type  <Return>   to continue : ")

x = list(a = 1:10, beta = exp(-3:3),
logic = c(TRUE,FALSE,FALSE,TRUE))
print(lapply(x,mean))
print(sapply(x,mean))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 1.4.5

attach(faithful) #resident dataset
plot(faithful)

S=readline(prompt="Type  <Return>   to continue : ")

plot(as.vector(time(mdeaths)),as.vector(mdeaths),cex=.6,pch=19,xlab="",ylab="Monthly deaths from bronchitis")
lines(spline(mdeaths),lwd=2,col="red",lty=3)
ar=arima(mdeaths,order=c(1,0,0))$coef
lines(as.vector(time(mdeaths))[-1], ar[2]+ar[1]*(mdeaths[-length(mdeaths)]-ar[2]),col="blue",lwd=2,lty=2)
title("Splines versus AR(1) predictor")
legend(1974,2800,legend=c("spline","AR(1)"),col=c("red","blue"),lty=c(3,2),lwd=c(2,2),cex=.5)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 1.7

sqrnt=function(y) {
  x=y/2
  while (abs(x*x-y) > 1e-10) x=(x+y/x)/2
  x
  }

sqrnt(9)

sqrnt(2)
