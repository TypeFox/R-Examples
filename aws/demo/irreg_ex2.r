require(aws)
if(exists("X11")) X11(,12,4.5)
par(mfrow=c(1,3),mar=c(3,3,3,.25),mgp=c(2,1,0))
n <- as.integer(readline("Sample size: Enter for 1000, provide positive integer (>500) otherwise"))
if(is.na(as.numeric(n))) n <- 1000 else n <- max(500,as.integer(n))
set.seed(1)
x <- runif(n,-1,1)
y <- runif(n,-1,1)
fbi<-function(x,y,k,r){
z1<-sqrt(x^2+y^2)
theta<-asin(x/z1)
z<-sin(k*theta)
z[z1<r]<-sin(pi*z1/r)[z1<r]
sign((x+y)*(x-y))*sign(z)
}
k <- as.integer(readline("Select number of waves: Enter for 5, provide positive integer otherwise "))
if(is.na(as.numeric(k))) k <- 5 else k <- max(1,as.integer(k))
z0 <- fbi(x,y,k,.5)
sigma <- readline("Standard deviation of noise:\n Press 'Enter' for sigma=.25, otherwise provide value of sigma:")
if(is.na(as.numeric(sigma))) sigma <- .25 else sigma <- as.numeric(sigma)
if(sigma <= 0) sigma <- .25
z <- rnorm(z0,z0,sigma)
nbins <- readline("Number of bins:\n Press 'Enter' for nbins=250, otherwise provide number of bins:")
if(is.na(as.numeric(nbins))) nbins <- 250 else nbins <- as.numeric(nbins)
if(nbins <= 10) nbins <- 125
hmax <- readline("Maximal bandwidth:\n Press 'Enter' for hmax=25, otherwise provide value of hmax:")
if(is.na(as.numeric(hmax))) hmax <- 25 else hmax <- as.numeric(hmax)
if(hmax <= 1) hmax <- 25
memory <- if(readline("Memory control(N/Y) :") %in% c("y","Y")) TRUE else FALSE
if(k<25) {
yhat <- aws.irreg(z,cbind(x,y),hmax=hmax,graph=TRUE,memory=memory,nbins=nbins,henv=15)
} else {
yhat <- aws.irreg(z,cbind(x,y),hmax=hmax,graph=TRUE,memory=memory,nbins=nbins,sigma2=sigma^2,henv=15)
}
readline("Press ENTER to show results")
oldpar <- par(mfrow=c(1,3),mar=c(3,3,3,.25),mgp=c(2,1,0))
im0 <- outer(seq(-1,1,length=nbins),seq(-1,1,length=nbins),"fbi",k,.5)
image(im0,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Artificial image (k=",k,")"))
image(awsdata(yhat,"y"),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Noisy image (sigma=",signif(sigma,3),")"))
image(awsdata(yhat,"est"),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  hmax=",signif(yhat@hmax,3)))
readline("Press ENTER to see next plot")
image(awsdata(yhat,"sd"),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Standard deviation of estimates (min:",signif(min(awsdata(yhat,"sd"),na.rm=TRUE),3)," max:",signif(max(awsdata(yhat,"sd"),na.rm=TRUE),3),")"))
image(awsdata(yhat,"ni"),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Number of observations per bin (range = 0 - ",max(awsdata(yhat,"ni")),")"))
image(awsdata(yhat,"mask"),zlim=c(0,1),col=gray((0:255)/255),xaxt="n",yaxt="n")
title("Mask of evaluated bins")

par(oldpar)
if(! readline("keep files and device (N/Y) :") %in% c("y","Y")){ 
rm(fbi,x,im0,y,z,sigma,hmax,yhat,memory)
dev.off()
}
