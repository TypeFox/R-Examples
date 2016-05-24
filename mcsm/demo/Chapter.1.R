# Chapter 1 R commands

# Section 1.3.2

x=matrix(1:4,ncol=3)
print(x[x>5])
print(x[1.])

S=readline(prompt="Type  <Return>   to continue : ")

# Section 1.3.3

x = list(a = 1:10, beta = exp(-3:3),
logic = c(TRUE,FALSE,FALSE,TRUE))
print(lapply(x,mean))
print(sapply(x,mean))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 1.4

x=rnorm(20)
y=3*x+5+rnorm(20,sd=0.3)
reslm=lm(y~x)
print(summary(reslm))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 1.5

b=matrix(1:9,ncol=3)
print(var(b))
print(sd(b)^2)
x=rnorm(25) #produces a N(0,1) sample of size 25
out=t.test(x)
print(names(out))

S=readline(prompt="Type  <Return>   to continue : ")

attach(faithful) #resident dataset
print(cor.test(faithful[,1],faithful[,2]))
print(ks.test(jitter(faithful[,1]),pnorm))
print(shapiro.test(faithful[,2]))
print(wilcox.test(faithful[,1]))

#S=readline(prompt="Type  <Return>   to continue : ")

x=seq(-3,3,le=5) # equidispersed regressor
y=2+4*x+rnorm(5) # simulated variable
print(lm(y~x))
print(summary(lm(y~x)))
out=lm(y~x)
print(sqrt(sum(out$res^2)/out$df))

#S=readline(prompt="Type  <Return>   to continue : ")

print(summary(lm(weight ~ feed, data = chickwts)))
print(anova(lm(weight ~ feed, data = chickwts)))

S=readline(prompt="Type  <Return>   to continue : ")

library(MASS)
print(glm(formula = type ~ bmi + age, family = "binomial", data = Pima.tr))
print(arima(diff(EuStockMarkets[,1]),order=c(0,0,5)))
print(acf(ldeaths, plot=F))
print(acf(ldeaths))
print(acf(ldeaths,type="partial"))

S=readline(prompt="Type  <Return>   to continue : ")

# Section 1.6

plot(faithful)

jpeg(file="faith.jpg")
par(mfrow=c(1,2),mar=c(4,2,2,1))
hist(faithful[,1],nclass=21,col="grey",main="",xlab=names(faithful)[1])
hist(faithful[,2],nclass=21,col="wheat",main="",xlab=names(faithful)[2])
dev.off()

S=readline(prompt="Type  <Return>   to continue : ")

plot(as.vector(time(mdeaths)),as.vector(mdeaths),cex=.6,pch=19,xlab="",ylab="Monthly deaths from bronchitis in the UK")
lines(spline(mdeaths),lwd=2,col="chocolate",lty=3)
ar=arima(mdeaths,order=c(1,0,0))$coef
lines(as.vector(time(mdeaths))[-1], ar[2]+ar[1]*(mdeaths[-length(mdeaths)]-ar[2]),col="grey",lwd=2,lty=2)
title("Splines versus AR(1) predictive")
ar=arima(mdeaths,order=c(1,0,0),seasonal=list(order=c(1,0,0),period=12))$coef
lines(as.vector(time(mdeaths))[-(1:13)],ar[3]+ar[1]*(mdeaths[-c(1:12,72)]-ar[3])+ar[2]*(mdeaths[-(60:72)]-ar[3]),
lwd=2,col="steelblue",lty=2)
title("\n\nand SAR(1,12) predictive")
legend(1974,2800,legend=c("spline","AR(1)","SAR(1,12)"),col=c("chocolate","grey","steelblue"),
lty=c(3,2,2),lwd=rep(2,3),cex=.5)

S=readline(prompt="Type  <Return>   to continue : ")

# Section 1.7

sqrnt=function(y) {

  x=y/2
  while (abs(x*x-y) > 1e-10) x=(x+y/x)/2
  x
  }

sqrnt(9)

sqrnt(2)
