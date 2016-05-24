samp.dist.mech <- function(rep, int=.05){
old.par <- par(no.readonly = TRUE)

g1 <- NULL; rm(g1)
g1 <- read.pnm(system.file("pictures/goat1.pgm", package="asbio"))

ymx <- 0.5; ymn <- 3
xback <- rnorm(1000,90.5,15); xmx <- max(xback); xmn <- min(xback)

retry <- function(){ 
xback <- rnorm(1000,90.5,15); xmx <- max(xback); xmn <- min(xback)
}

if(xmx > 110|xmn < 75)retry() 
if(xmx > 110|xmn < 75)retry()
                             

m <- (ymx - ymn)/(xmx - xmn)
b <- -m*xmx + ymx
wt1 <- matrix(nrow=10,ncol=rep,data=sample(xback, rep*10, replace=FALSE))
all.wt <- apply(wt1,2,mean)


sub <- function(wt){
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11),5,5, byrow= TRUE))
fit <- b + m * wt[1]; par(mar=rep(fit, 4))
plot(g1); mtext(side=1,line= .8,bquote(paste(.(round(wt[1],1)),"kg")))
fit <- b + m * wt[2]; par(mar=rep(fit, 4))
plot(g1); mtext(side=1,line= .8,bquote(paste(.(round(wt[2],1)),"kg")))
fit <- b + m * wt[3]; par(mar=rep(fit, 4))
plot(g1); mtext(side=1,line= .8,bquote(paste(.(round(wt[3],1)),"kg")))
fit <- b + m * wt[4]; par(mar=rep(fit, 4))
plot(g1); mtext(side=1,line= .8,bquote(paste(.(round(wt[4],1)),"kg")))
fit <- b + m * wt[5]; par(mar=rep(fit, 4))
plot(g1); mtext(side=1,line= .8,bquote(paste(.(round(wt[5],1)),"kg")))
fit <- b + m * wt[6]; par(mar=rep(fit, 4))
plot(g1); mtext(side=1,line= .8,bquote(paste(.(round(wt[6],1)),"kg")))
fit <- b + m * wt[7]; par(mar=rep(fit, 4))
plot(g1); mtext(side=1,line= .8,bquote(paste(.(round(wt[7],1)),"kg")))
fit <- b + m * wt[8]; par(mar=rep(fit, 4))
plot(g1); mtext(side=1,line= .8,bquote(paste(.(round(wt[8],1)),"kg")))
fit <- b + m * wt[9]; par(mar=rep(fit, 4))
plot(g1); mtext(side=1,line= .8,bquote(paste(.(round(wt[9],1)),"kg")))
fit <- b + m * wt[10]; par(mar=rep(fit, 4))
plot(g1); mtext(side=1,line= .8,bquote(paste(.(round(wt[10],1)),"kg")))
}



for(i in 1:rep){
dev.hold()
sub(wt1[,i])
mtext(side = 3, outer=TRUE,"", line = -26, cex = 1.3)
Sys.sleep(int*10)
mtext(side = 3, outer=TRUE,bquote(paste(italic(bar(x))[.(i)], " = ", .(round(all.wt[i],1)))), line = -24, cex = 1.3)

par(mar=c(6,5.5,5,1))
hist(all.wt[1:i], main = "", breaks = c(75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110), cex.lab=1.8, cex.axis=1.5, xlab = "", ylab = "", ylim=c(0,32))
arrows(all.wt[i], 35, all.wt[i], 28, length = 0.25, lwd = 10, col="black", angle= 30)
mtext(side=1, expression(italic(bar(x))), cex=1.3, line = 4)
mtext(side=2, "Cumulative freq.", cex = 1.3, line = 4)
dev.flush()
Sys.sleep(int)
}
on.exit(par(old.par))
}

