require(aws)
if(exists("X11")) X11(,12,4.5)
f1 <- function(x){
#  
#    Blocks data (Example 6 from Fan & Gijbels (1996)
#
     xj <- c(10,13,15,23,25,40,44,65,76,78,81)/100
     hj <- c(40,-50,30,-40,50,-42,21,43,-31,21,-42)*.37
     Kern <- function(x) (1-sign(x))/2
     apply(Kern(outer(xj,x,"-"))*hj,2,sum)
     }
f2 <- function(x) 25*sin(2*pi*1.2/(x+.2))

x <- seq(0,1,length=2048)
example <- as.integer(readline("Select example: 1 for Blocks data (default),  2 for Smooth example"))
if(!(example %in% (1:2))) example <- 1 
u <- switch(example,f1(x),f2(x))
plot(x,u,type="l",col=2,lwd=2)
title("Regression function")
sigma <- readline("Standard deviation of noise:\n Press 'Enter' for sigma=3, otherwise provide value of sigma:")
if(is.na(as.numeric(sigma))) sigma <- 3 else sigma <- as.numeric(sigma)
if(sigma <= 0) sigma <- 3
y <- rnorm(u,u,sigma)
plot(x,u,type="l",col=2,ylim=range(u,y),lwd=2)
title("Regression function and data")
points(x,y)
hmax <- readline("Maximal bandwidth:\n Press 'Enter' for hmax=250, otherwise provide value of hmax:")
if(is.na(as.numeric(hmax))) hmax <- 250 else hmax <- as.numeric(hmax)
if(hmax <= 1) hmax <- 250
cat("Run aws \n")
yhat <- aws(y,hmax=hmax,graph=TRUE)
readline("Press ENTER to show results")
plot(x,u,type="l",col=3,ylim=range(u,y),lwd=2,lty=2)
points(x,y)
lines(x,awsdata(yhat,"est"),col=2,lwd=2)
title("Data, fitted (red) and true (green) values")
if(! readline("keep files and device (N/Y) :") %in% c("y","Y")){ 
rm(f1,f2,x,u,y,sigma,hmax,yhat)
dev.off()
}
