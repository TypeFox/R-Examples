require(aws)
if(exists("X11")) X11(,12,4.5)
par(mfrow=c(1,3),mar=c(3,3,3,.25),mgp=c(2,1,0))
x<-seq(-1,1,length=256)
fbi<-function(x,y,k,r){
z1<-sqrt(x^2+y^2)
theta<-asin(x/z1)
z<-sin(k*theta)
z[z1<r]<-sin(pi*z1/r)[z1<r]
z<-sign((x+y)*(x-y))*z
(z-min(z))/(max(z)-min(z))
}
k <- as.integer(readline("Select number of waves: Enter for 7, provide positive integer otherwise "))
if(is.na(as.numeric(k))) k <- 7 else k <- max(1,as.integer(k))
im0<-outer(x,x,"fbi",k,.5)
image(im0,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Artificial image (k=",k,")"))
sigma <- readline("Standard deviation of noise:\n Press 'Enter' for sigma=.25, otherwise provide value of sigma:")
if(is.na(as.numeric(sigma))) sigma <- .25 else sigma <- as.numeric(sigma)
if(sigma <= 0) sigma <- .25
y <- im0+rnorm(im0,0,sigma)
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Noisy image (sigma=",signif(sigma,3),")"))
degree <- readline("Degree of polynomial model:\n Press 'Enter' for degree =2, otherwise provide degree:")
if(is.na(as.numeric(degree))) degree <- 2 else degree <- as.numeric(degree)
if(!(degree %in% 0:2)) degree <- 2
hmax <- readline("Maximal bandwidth:\n Press 'Enter' for hmax=15, otherwise provide value of hmax:")
if(is.na(as.numeric(hmax))) hmax <- 15 else hmax <- as.numeric(hmax)
if(hmax <= 1) hmax <- 15
memory <- if(readline("Memory control(N/Y) :") %in% c("y","Y")) TRUE else FALSE
risk <- readline("Report risks (N/Y):")
if(risk %in% c("y","Y")) u <-im0 else u <- NULL
cat("Run aws \n")
if(k<25) {
yhat <- lpaws(y,degree=degree,hmax=hmax,graph=TRUE,memory=memory,u=u)
} else {
yhat <- lpaws(y,degree=degree,hmax=hmax,graph=TRUE,memory=memory,sigma2=sigma^2,u=u)
}
readline("Press ENTER to show results")
oldpar <- par(mfrow=c(1,3),mar=c(3,3,3,.25),mgp=c(2,1,0))
image(im0,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Artificial image (k=",k,")"))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Noisy image (sigma=",signif(sigma,3),")"))
image(awsdata(yhat,"est"),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction (degree=",yhat@degree," hmax=",signif(yhat@hmax,3),")"))
readline("Enter for next plot:")
image(awsdata(yhat,"sd"),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Standard deviation of estimates (min:",signif(min(awsdata(yhat,"sd")),3)," max:",signif(max(awsdata(yhat,"sd")),3),")"))
if(degree>0) image(awsdata(yhat,"theta")[,,2],col=gray((0:255)/255),
                   xaxt="n",yaxt="n",main="first derivative")
if(degree>1) image(awsdata(yhat,"theta")[,,3],col=gray((0:255)/255),
                   xaxt="n",yaxt="n",main="second derivative")
par(oldpar)
if(! readline("keep files and device (N/Y) :") %in% c("y","Y")){ 
rm(fbi,x,im0,y,sigma,hmax,yhat,degree,u,risk,memory)
dev.off()
}
