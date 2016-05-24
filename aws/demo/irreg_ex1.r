require(aws)
if(exists("X11")) X11(,10,10)
n <- as.integer(readline("Sample size: Enter for 1000, provide positive integer (>250) otherwise"))
if(is.na(as.numeric(n))) n <- 1000 else n <- max(250,as.integer(n))
set.seed(1)
x <- runif(n,0,1)
f1 <- function(x){
#  
#    Blocks data (Example 6 from Fan & Gijbels (1996)
#
     xj <- c(10,13,15,23,25,40,44,65,76,78,81)/100
     hj <- c(40,-50,30,-40,50,-42,21,43,-31,21,-42)*.37
     Kern <- function(x) (1-sign(x))/2
     apply(Kern(outer(xj,x,"-"))*hj,2,sum)
     }
y0 <- f1(x)
sigma <- readline("Standard deviation of noise:\n Press 'Enter' for sigma=3, otherwise provide value of sigma:")
if(is.na(as.numeric(sigma))) sigma <- 3 else sigma <- as.numeric(sigma)
if(sigma <= 0) sigma <- 3
y <- rnorm(y0,y0,sigma)
nbins <- readline("Number of bins:\n Press 'Enter' for nbins=2000, otherwise provide number of bins:")
if(is.na(as.numeric(nbins))) nbins <- 2000 else nbins <- as.numeric(nbins)
if(nbins <= 250) nbins <- 250
hmax <- readline("Maximal bandwidth:\n Press 'Enter' for hmax=500, otherwise provide value of hmax:")
if(is.na(as.numeric(hmax))) hmax <- 500 else hmax <- as.numeric(hmax)
if(hmax <= 1) hmax <- 50
memory <- if(readline("Memory control(N/Y) :") %in% c("y","Y")) TRUE else FALSE
yhat <- aws.irreg(y,x,hmax=hmax,graph=TRUE,memory=memory,nbins=nbins,henv=1+3*nbins/n)
readline("Press ENTER to show results")
oldpar <- par(mfrow=c(1,1),mar=c(3,3,3,.25),mgp=c(2,1,0))
plot(x,y)
lines(seq(yhat@xmin,yhat@xmax,length=yhat@dy),awsdata(yhat,"est"),col=2,lwd=3)
lines(seq(yhat@xmin,yhat@xmax,length=yhat@dy),f1(seq(yhat@xmin,yhat@xmax,length=yhat@dy)),lty=2,lwd=3,col=3)
title("Observed data, reconstruction (red) and true function (green)")
par(oldpar)
if(! readline("keep files and device (N/Y) :") %in% c("y","Y")){ 
rm(f1,x,y,sigma,hmax,yhat,memory)
dev.off()
}
