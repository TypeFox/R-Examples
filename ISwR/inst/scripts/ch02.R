load("ch1.RData")  
.foo <- dev.copy2eps
rm(dev.copy2eps)
ls()
dev.copy2eps <- .foo
rm(height, weight)
sink("myfile")
ls()
sink()
attach(thuesen)
blood.glucose
search()
detach()
search()
thue2 <- subset(thuesen,blood.glucose<7)
thue2
thue3 <- transform(thuesen,log.gluc=log(blood.glucose))
thue3
thue4 <- within(thuesen,{
   log.gluc <- log(blood.glucose)
   m <- mean(log.gluc)
   centered.log.gluc <- log.gluc - m
   rm(m)
})
thue4
d <- par(mar=c(5,4,4,2)+.1)
x <- runif(50,0,2)
y <- runif(50,0,2)
plot(x, y, main="Main title", sub="subtitle",
     xlab="x-label", ylab="y-label")
text(0.6,0.6,"text at (0.6,0.6)")
abline(h=.6,v=.6)
for (side in 1:4) mtext(-1:4,side=side,at=.7,line=-1:4)
mtext(paste("side",1:4), side=1:4, line=-1,font=2)
if (.make.epsf) dev.copy2eps(file="layout.ps")
par(d)
plot(x, y, type="n", xlab="", ylab="", axes=F)
points(x,y)
axis(1)
axis(2,at=seq(0.2,1.8,0.2))
box()
title(main="Main title", sub="subtitle",
    xlab="x-label", ylab="y-label")
set.seed(1234) #make it happen....
x <- rnorm(100)
hist(x,freq=F)
curve(dnorm(x),add=T)  
h <- hist(x, plot=F)
ylim <- range(0, h$density, dnorm(0))
hist(x, freq=F, ylim=ylim)
curve(dnorm(x), add=T)  
if (.make.epsf) dev.copy2eps(file="hist+norm.ps")
hist.with.normal <- function(x, xlab=deparse(substitute(x)),...)
{
    h <- hist(x, plot=F, ...)
    s <- sd(x)
    m <- mean(x)
    ylim <- range(0,h$density,dnorm(0,sd=s))
    hist(x, freq=F, ylim=ylim, xlab=xlab, ...)
    curve(dnorm(x,m,s), add=T)
}
hist.with.normal(rnorm(200))
y <- 12345
x <- y/2 
while (abs(x*x-y) > 1e-10) x <- (x + y/x)/2
x
x^2  
x <- y/2 
repeat{ 
    x <- (x + y/x)/2
    if (abs(x*x-y) < 1e-10) break
}
x
x <- seq(0, 1,.05)
plot(x, x, ylab="y", type="l")
for ( j in 2:8 ) lines(x, x^j)
t.test(bmi, mu=22.5)$p.value
print
length(methods("print")) # quoted in text
thuesen2 <- read.table(
   system.file("rawdata","thuesen.txt",package="ISwR"), header=T)
thuesen2
levels(secretin$time)
system.file("rawdata", "thuesen.txt", package="ISwR")
rm(list=ls())
while(search()[2] != "package:ISwR") detach()
