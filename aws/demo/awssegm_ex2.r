require(aws)
if(exists("X11")) X11(,10,10)
y <- matrix(rnorm(512^2),512,512)

y1 <- numeric(1024)
y1[-21+64:65] <- 1
y1[-21+127:130] <- sqrt(.5)
y1[-21+189:196] <- .5
y1[-21+249:264] <- sqrt(.125)
y1[-21+305:336] <- .25
y1[320+64:65] <-   sqrt(2)
y1[320+127:130] <- 1
y1[320+189:196] <- sqrt(.5)
y1[320+249:264] <- .5
y1[320+305:336] <- sqrt(.25)
y1[661+64:65] <- 2
y1[661+127:130] <- sqrt(2) 
y1[661+189:196] <- 1
y1[661+249:264] <- sqrt(.5)
y1[661+305:336] <- .5
y1 <- outer(y1,y1,"*")

y1 <- numeric(256)
y1[10:11] <- 1
y1[32:35] <- sqrt(.5)
y1[56:63] <- .5
y1[84:99] <- sqrt(.125)
y1[120:151] <- .25
y1[172:235] <- sqrt(.03125)
y1 <- outer(y1,y1,"*")
y1 <- cbind(y1,sqrt(2)*y1)
y1 <- rbind(y1,2*y1)

delta <- as.numeric(readline("Width of central interval (default 0)"))
if(is.na(delta)||delta <= 0) delta <- 0 else delta <- delta/2

size <- as.numeric(readline("Signal size (default 5+intervalwidth/2)"))
if(is.na(size)||size <= 0) size <- 5+delta

par(mfrow=c(2,2),mar=c(3,3,3,1),mgp=c(2,1,0))
image(size*y1,col=grey((0:255)/255))
title("Original image")
image(y+size*y1,col=grey((0:255)/255))
title("Noisy image")
truesegm <- (size*y1>delta)-(size*y1<delta)
image(truesegm)
title("True segmentation")
readline("Enter for next step")

z<-aws.segment(y+size*y1,level=0,delta=delta,hmax=40,graph=TRUE)
dim(z@theta) <- dim(z@theta)[1:2]
ergs<-c(sum(z@segment>0&size*y1>delta),sum(z@segment>0&size*y1<=delta),sum(z@segment<=0&size*y1>delta),sum(z@segment<=0&size*y1<=delta))
names(ergs) <-c("Correct alternative","False positive","False negative","Correct hypothesis")
print(ergs)
readline("Enter for next step")

par(mfrow=c(2,2),mar=c(3,3,3,1),mgp=c(2,1,0))
image(size*y1,col=grey((0:255)/255))
title("Original image")
image(y+size*y1,col=grey((0:255)/255))
title("Noisy image")
image(truesegm) 
title("True segmentation")
image(z@segment) 
title("Segmentation results")

if(! readline("keep files and device (N/Y) :") %in% c("y","Y")){ 
rm(y,y1,size,delta,z,ergs,truesegm)
dev.off()
}
