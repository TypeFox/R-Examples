library(aws)
y <- array(rnorm(32^2*20),c(32,32,20))

y1 <- array(0,c(32,32,20))
y1[5:6,5:6,3:4] <- 1
y1[11:14,5:6,3:4] <- sqrt(.5)
y1[21:28,5:6,3:4] <- .5
y1[5:6,11:14,3:4] <- sqrt(.5)
y1[11:14,11:14,3:4] <- .5
y1[21:28,11:14,3:4] <- sqrt(.125)
y1[5:6,21:28,3:4] <- .5
y1[11:14,21:28,3:4] <- sqrt(.125)
y1[21:28,21:28,3:4] <- .25

y1[5:6,5:6,8:10] <- 1
y1[11:14,5:6,8:10] <- sqrt(.5)
y1[21:28,5:6,8:10] <- .5
y1[5:6,11:14,8:10] <- sqrt(.5)
y1[11:14,11:14,8:10] <- .5
y1[21:28,11:14,8:10] <- sqrt(.125)
y1[5:6,21:28,8:10] <- .5
y1[11:14,21:28,8:10] <- sqrt(.125)
y1[21:28,21:28,8:10] <- .25

y1[5:6,5:6,15:18] <- 1
y1[11:14,5:6,15:18] <- sqrt(.5)
y1[21:28,5:6,15:18] <- .5
y1[5:6,11:14,15:18] <- sqrt(.5)
y1[11:14,11:14,15:18] <- .5
y1[21:28,11:14,15:18] <- sqrt(.125)
y1[5:6,21:28,15:18] <- .5
y1[11:14,21:28,15:18] <- sqrt(.125)
y1[21:28,21:28,15:18] <- .25


delta <- as.numeric(readline("Width of central interval (default 4)"))
if(is.na(delta)||delta < 0) delta <- 4 
delta <- delta/2

size <- as.numeric(readline("Signal size (default 5+intervalwidth/2)"))
if(is.na(size)||size <= 0) size <- 5+delta

if(exists("X11")) X11(,15,6.5)

par(mfrow=c(4,10),mar=c(2,2,2,1),mgp=c(2,1,0))
for(i in 1:20) image(y[,,i]+size*y1[,,i],col=grey((0:255)/255),main=paste("Noisy slice",i))
for(i in 1:20) image(size*y1[,,i]>delta,col=grey((0:255)/255),main=paste("True segm.",i))
readline("Enter for next step")


z<-aws.segment(y+size*y1,level=0,delta=delta,hmax=10,graph=FALSE,ladjust=1,ext=0)
dim(z@theta) <- dim(z@theta)[1:3]
ergs<-c(sum(z@segment>0&size*y1>delta),sum(z@segment>0&size*y1<=delta),sum(z@segment<=0&size*y1>delta),sum(z@segment<=0&size*y1<=delta))
names(ergs) <-c("Correct alternative","False positive","False negative","Correct hypothesis")
print(ergs)
readline("Enter for next step")

if(exists("X11")) X11(,15,6.5)

truesegm <- (size*y1>delta)-(size*y1<delta)
par(mfrow=c(4,10),mar=c(2,2,2,1),mgp=c(2,1,0))
for(i in 1:20) image(z@segment[,,i],col=grey((0:255)/255),main=paste("Segm.res.",i))
for(i in 1:20) image(truesegm[,,i],col=grey((0:255)/255),main=paste("True segm.",i))

if(! readline("keep files and device (N/Y) :") %in% c("y","Y")){ 
rm(y,y1,size,delta,z,ergs,truesegm)
dev.off();dev.off()
}
