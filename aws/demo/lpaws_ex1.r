require(aws)
if(exists("X11")) X11(,10,5)
f1 <- function(x){
   12*sin(2*pi*x)+3*sign(x>.4)-3*sign(x>.7)
     }
f2 <- function(x) 10*sin(2*pi*1.2/(x+.2))

x <- seq(0,1,length=2048)
example <- as.integer(readline("Select example: 1 for discontinuous (default),  2 for smooth example"))
if(!(example %in% (1:2))) example <- 1 
u <- switch(example,f1(x),f2(x))
plot(x,u,type="l",col=2,lwd=2)
title("Regression function")
sigma <- readline("Standard deviation of noise:\n Press 'Enter' for sigma=1, otherwise provide value of sigma:")
if(is.na(as.numeric(sigma))) sigma <- 1 else sigma <- as.numeric(sigma)
if(sigma <= 0) sigma <- 1
y <- rnorm(u,u,sigma)
plot(x,u,type="l",col=2,ylim=range(u,y),lwd=2)
title("Regression function and data")
points(x,y)
degree <- readline("Degree of polynomial model:\n Press 'Enter' for degree =2, otherwise provide degree:")
if(is.na(as.numeric(degree))) degree <- 2 else degree <- as.numeric(degree)
if(!(degree %in% 0:2)) degree <- 2
hmax <- readline("Maximal bandwidth:\n Press 'Enter' for hmax=300, otherwise provide value of hmax:")
if(is.na(as.numeric(hmax))) hmax <- 250 else hmax <- as.numeric(hmax)
if(hmax <= 1) hmax <- 250
memory <- if(readline("Memory control(N/Y) :") %in% c("y","Y")) TRUE else FALSE
risk <- readline("Report risks (N/Y):")
if(risk %in% c("y","Y")) uu <-u else uu <- NULL
cat("Run aws \n")
yhat <- lpaws(y,degree=degree,hmax=hmax,graph=TRUE,memory=memory,u=uu)
readline("Press ENTER to show results")
plot(x,u,type="l",col=2,ylim=range(u,y,awsdata(yhat,"est")),lwd=2)
points(x,y)
lines(x,awsdata(yhat,"est"),col=2,lwd=2)
lines(x,u,col=3,lty=2,lwd=2)
title("Data, fitted (red) and true (green) values")
if(! readline("keep files and device (N/Y) :") %in% c("y","Y")){ 
rm(f1,f2,x,u,y,sigma,hmax,yhat,degree,uu,risk,memory)
dev.off()
}
