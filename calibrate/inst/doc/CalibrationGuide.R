### R code from vignette source 'CalibrationGuide.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: CalibrationGuide.Rnw:79-92
###################################################
prefig <- function(){
data(goblets)
X <- goblets
plot(X[,1],X[,2],pch=19,cex=0.5,xlab=expression(X[1]),
ylab=expression(X[2]),xlim=c(5,25),ylim=c(5,25),asp=1)
m <- apply(cbind(X[,1],X[,2]),2,mean)
textxy(X[,1],X[,2],1:25,m=m,cex=0.75)
points(m[1],m[2],col="red",pch=19,cex=0.5)
abline(h = m[2])
abline(v = m[1])
}
options("SweaveHooks"=list(aap=prefig))
options("width"=60)


###################################################
### code chunk number 2: noot
###################################################
library(calibrate)


###################################################
### code chunk number 3: CalibrationGuide.Rnw:113-115
###################################################
data(goblets)
X <- goblets


###################################################
### code chunk number 4: CalibrationGuide.Rnw:125-130
###################################################
plot(X[,1],X[,2],pch=19,cex=0.5,xlab=expression(X[1]),ylab=expression(X[2]),
     xlim=c(5,25),ylim=c(5,25),asp=1)
m <- apply(X[,1:2],2,mean)
textxy(X[,1],X[,2],1:25,m=m,cex=0.75)
origin(m)


###################################################
### code chunk number 5: CalibrationGuide.Rnw:137-143
###################################################
getOption("SweaveHooks")[["aap"]]()
Xc <- scale(X,center=TRUE,scale=FALSE)
b <- solve(t(Xc[,1:2])%*%Xc[,1:2])%*%t(Xc[,1:2])%*%Xc[,5]
print(b)
bscaled <- 20*b
arrows(m[1],m[2],m[1]+bscaled[1],m[2]+bscaled[2],col="blue",length=0.1)
arrows(m[1],m[2],m[1]-bscaled[1],m[2]-bscaled[2],length=0,lty="dashed",col="blue")


###################################################
### code chunk number 6: CalibrationGuide.Rnw:155-161
###################################################
getOption("SweaveHooks")[["aap"]]()
print(range(X[,5]))
yc <- scale(X[,5],scale=FALSE)
tm <- seq(2,10,by=1)
tmc <- tm - mean(X[,5])
Calibrate.X5<-calibrate(b,yc,tmc,Xc[,1:2],tmlab=tm,m=m,tl=0.3,axislab="X_5",
                        labpos=4,cex.axislab=1)


###################################################
### code chunk number 7: CalibrationGuide.Rnw:198-203
###################################################
getOption("SweaveHooks")[["aap"]]()
yc <- scale(X[,5],scale=FALSE)
tm <- seq(2,10,by=1)
tmc <- tm - mean(X[,5])
Calibrate.X5<-calibrate(b,yc,tmc,Xc[,1:2],tmlab=tm,m=m,tl=0.3,axislab="X_5",labpos=4,
                        cex.axislab=1,shiftdir="right")


###################################################
### code chunk number 8: CalibrationGuide.Rnw:220-229
###################################################
getOption("SweaveHooks")[["aap"]]()
opar <- par('xpd'=TRUE)
tm <- seq(3,8,by=1)
tmc <- (tm - mean(X[,5]))
Calibrate.rightmargin.X5 <- calibrate(c(0,1),yc,tmc,Xc[,1:2],tmlab=tm,m=m,
                                      axislab="X_5",tl=0.5,
                                      shiftvec=c(par('usr')[2]-mean(X[,1]),0),
                                      shiftfactor=1,where=2,
                                      laboffset=c(1.5,1.5),cex.axislab=1)
par(opar)


###################################################
### code chunk number 9: CalibrationGuide.Rnw:252-260
###################################################
getOption("SweaveHooks")[["aap"]]()
tm <- seq(2,10,by=1)
tmc <- (tm - mean(X[,5]))
Calibrate.X5 <- calibrate(b,yc,tmc,Xc[,1:2],tmlab=tm,m=m,axislab="X_5",tl=0.5,
                          dp=TRUE,labpos=4)
tm <- seq(2,10,by=0.1)
tmc <- (tm - mean(X[,5]))
Calibrate.X5 <- calibrate(b,yc,tmc,Xc[,1:2],tmlab=tm,m=m,tl=0.25,verb=FALSE,
                          lm=FALSE)


###################################################
### code chunk number 10: CalibrationGuide.Rnw:276-298
###################################################
getOption("SweaveHooks")[["aap"]]()
opar <- par('xpd'=TRUE)
tm <- seq(5,25,by=5)
tmc <- (tm - mean(X[,1]))
yc <- scale(X[,1],scale=FALSE)
Calibrate.X1 <- calibrate(c(1,0),yc,tmc,Xc[,1:2],tmlab=tm,m=m,tl=0.5,
                          axislab="X_1",cex.axislab=1,showlabel=FALSE,
shiftvec=c(0,-(m[2]-par("usr")[3])),shiftfactor=1,reverse=TRUE)
tm <- seq(5,25,by=1); tmc <- (tm - mean(X[,1]))
Calibrate.X1 <- calibrate(c(1,0),yc,tmc,Xc[,1:2],tmlab=tm,m=m,tl=0.25,
                          axislab="X_1",cex.axislab=1,showlabel=FALSE,
shiftvec=c(0,-(m[2]-par("usr")[3])),shiftfactor=1,reverse=TRUE,
                          verb=FALSE,lm=FALSE)
yc <- scale(X[,2],scale=TRUE)
tm <- seq(-3,1,by=1)
Calibrate.X2 <- calibrate(c(0,1),yc,tm,Xc[,1:2],tmlab=tm,m=m,tl=0.6,
                          axislab="X_2",cex.axislab=1,showlabel=FALSE,
shiftvec=c(-(mean(X[,1])-par('usr')[1]),0),shiftfactor=1,verb=TRUE,lm=TRUE)
tm <- seq(-3,1.5,by=0.1)
Calibrate.X2 <- calibrate(c(0,1),yc,tm,Xc[,1:2],tmlab=tm,m=m,tl=0.3,
                          axislab="X_2",cex.axislab=1,showlabel=FALSE,
shiftvec=c(-(mean(X[,1])-par('usr')[1]),0),shiftfactor=1,verb=FALSE,lm=FALSE)
par(opar)


###################################################
### code chunk number 11: CalibrationGuide.Rnw:327-346
###################################################
# PCA and Biplot construction
pca.results <- princomp(X,cor=FALSE)
Fp <- pca.results$scores
Gs <- pca.results$loadings
plot(Fp[,1],Fp[,2],pch=16,asp=1,xlab="PC 1",ylab="PC 2",cex=0.5)
textxy(Fp[,1],Fp[,2],rownames(X),cex=0.75)
arrows(0,0,15*Gs[,1],15*Gs[,2],length=0.1)
textxy(15*Gs[,1],15*Gs[,2],colnames(X),cex=0.75)
# Calibration of X_3
ticklab  <- seq(5,30,by=5)
ticklabc <- ticklab-mean(X[,3])
yc <- (X[,3]-mean(X[,3])) 
g <- Gs[3,1:2]                                  
Calibrate.X3 <- calibrate(g,yc,ticklabc,Fp[,1:2],ticklab,tl=0.5,
                          axislab="X3",cex.axislab=0.75,where=1,labpos=4)
ticklab  <- seq(5,30,by=1)
ticklabc <- ticklab-mean(X[,3])
Calibrate.X3.fine <- calibrate(g,yc,ticklabc,Fp[,1:2],ticklab,lm=FALSE,tl=0.25,
                               verb=FALSE,cex.axislab=0.75,where=1,labpos=4)


###################################################
### code chunk number 12: CalibrationGuide.Rnw:353-376
###################################################
# PCA and Biplot construction
pca.results <- princomp(X,cor=TRUE)
Fp <- pca.results$scores
Ds <- diag(pca.results$sdev)
Fs <- Fp%*%solve(Ds)
Gs <- pca.results$loadings
Gp <- Gs%*%Ds
#plot(Fs[,1],Fs[,2],pch=16,asp=1,xlab="PC 1",ylab="PC 2",cex=0.5)
#textxy(Fs[,1],Fs[,2],rownames(X))
plot(Gp[,1],Gp[,2],pch=16,cex=0.5,xlim=c(-1,1),ylim=c(-1,1),asp=1,
     xlab="1st principal axis",ylab="2nd principal axis")
arrows(0,0,Gp[,1],Gp[,2],length=0.1)
textxy(Gp[,1],Gp[,2],colnames(X),cex=0.75)
ticklab <- c(seq(-1,-0.2,by=0.2),seq(0.2,1.0,by=0.2))
R <- cor(X)
y <- R[,5]
g <- Gp[5,1:2]                                        
Calibrate.X5 <- calibrate(g,y,ticklab,Gp[,1:2],ticklab,lm=TRUE,tl=0.05,dp=TRUE,
                      labpos=2,cex.axislab=0.75,axislab="X_5")
ticklab <- seq(-1,1,by=0.1)
Calibrate.X5 <- calibrate(g,y,ticklab,Gp[,1:2],ticklab,lm=FALSE,tl=0.05,verb=FALSE)
ticklab <- seq(-1,1,by=0.01)
Calibrate.X5 <- calibrate(g,y,ticklab,Gp[,1:2],ticklab,lm=FALSE,tl=0.025,verb=FALSE)


###################################################
### code chunk number 13: CalibrationGuide.Rnw:389-391
###################################################
print(R)
print(cbind(R[,5],Calibrate.X5$yt))


###################################################
### code chunk number 14: CalibrationGuide.Rnw:430-472
###################################################
library(MASS)
data(calves)
ca.results <- corresp(calves,nf=2)
Fs <- ca.results$rscore
Gs <- ca.results$cscore
Ds <- diag(ca.results$cor)
Fp <- Fs%*%Ds
Gp <- Gs%*%Ds
plot(Gs[,1],Gs[,2],pch=16,asp=1,cex=0.5,xlab="1st principal axis",
     ylab="2nd principal axis")
textxy(Gs[,1],Gs[,2],colnames(calves),cex=0.75)
points(Fp[,1],Fp[,2],pch=16,cex=0.5)
textxy(Fp[,1],Fp[,2],rownames(calves),cex=0.75)
origin()
arrows(0,0,Gs[,1],Gs[,2])

P <- as.matrix(calves/sum(calves))
r <- apply(P,1,sum)
k <- apply(P,2,sum)
Dc <- diag(k)
Dr <- diag(r)

RP <- solve(Dr)%*%P
print(RP)
CRP <- RP - ones(nrow(RP), 1) %*% t(k)
TCRP <- CRP%*%solve(Dc)

y <- TCRP[,3]
g <- Gs[3,1:2]

ticklab  <- c(0,seq(0,1,by=0.2))                         
ticklabs <- (ticklab - k[3])/k[3]
Calibrate.AI <- calibrate(g,y,ticklabs,Fp[,1:2],ticklab,lm=TRUE,tl=0.10,
                          weights=Dr,axislab="AI",labpos=4,dp=TRUE)
ticklab  <- c(0,seq(0,1,by=0.1))                         
ticklabs <- (ticklab - k[3])/k[3]
Calibrate.AI <- calibrate(g,y,ticklabs,Fp[,1:2],ticklab,lm=FALSE,tl=0.10,
                          weights=Dr,verb=FALSE)
ticklab  <- c(0,seq(0,1,by=0.01))                         
ticklabs <- (ticklab - k[3])/k[3]
Calibrate.AI <- calibrate(g,y,ticklabs,Fp[,1:2],ticklab,lm=FALSE,tl=0.05,
                          weights=Dr,verb=FALSE)


###################################################
### code chunk number 15: CalibrationGuide.Rnw:490-521
###################################################
data(heads)
X  <- cbind(heads$X1,heads$X2)
Y  <- cbind(heads$Y1,heads$Y2)

Rxy<- cor(X,Y)
Ryx<- t(Rxy)
Rxx<- cor(X)
Ryy<- cor(Y)

cca.results <-canocor(X,Y)

plot(cca.results$Gs[,1],cca.results$Gs[,2],pch=16,asp=1,xlim=c(-1,1),ylim=c(-1,1),
     xlab=expression(V[1]),ylab=expression(V[2]))
arrows(0,0,cca.results$Fp[,1],cca.results$Fp[,2],length=0.1)
arrows(0,0,cca.results$Gs[,1],cca.results$Gs[,2],length=0.1)

textxy(cca.results$Fp[1,1],cca.results$Fp[1,2],expression(X[1]),cex=0.75)
textxy(cca.results$Fp[2,1],cca.results$Fp[2,2],expression(X[2]),cex=0.75)

textxy(cca.results$Gs[1,1],cca.results$Gs[1,2],expression(Y[1]),cex=0.75)
textxy(cca.results$Gs[2,1],cca.results$Gs[2,2],expression(Y[2]),cex=0.75)

circle(1)
ticklab  <- seq(-1,1,by=0.2)  

y <- Rxy[,2]
g <- cca.results$Gs[2,1:2]                        

Cal.Cor.Y2 <- calibrate(g,y,ticklab,cca.results$Fp[,1:2],ticklab,lm=TRUE,tl=0.05,
                        dp=TRUE,reverse=TRUE,weights=solve(Rxx),
axislab="Y_2",cex.axislab=0.75,showlabel=FALSE)


###################################################
### code chunk number 16: CalibrationGuide.Rnw:524-557
###################################################

plot(cca.results$Gs[,1],cca.results$Gs[,2],pch=16,asp=1,xlim=c(-2,2),ylim=c(-2,2),
     xlab=expression(V[1]),ylab=expression(V[2]))
#arrows(0,0,cca.results$Fp[,1],cca.results$Fp[,2],length=0.1)
#arrows(0,0,cca.results$Gs[,1],cca.results$Gs[,2],length=0.1)

textxy(cca.results$Fp[1,1],cca.results$Fp[1,2],expression(X[1]))
textxy(cca.results$Fp[2,1],cca.results$Fp[2,2],expression(X[2]))

textxy(cca.results$Gs[1,1],cca.results$Gs[1,2],expression(Y[1]))
textxy(cca.results$Gs[2,1],cca.results$Gs[2,2],expression(Y[2]))

points(cca.results$V[,1],cca.results$V[,2],pch=16,cex=0.5)
textxy(cca.results$V[,1],cca.results$V[,2],1:nrow(X),cex=0.75)


ticklab  <- seq(135,160,by=5)                           
ticklabc <- ticklab-mean(Y[,2])
ticklabs <- (ticklab-mean(Y[,2]))/sqrt(var(Y[,2]))

y <- (Y[,2]-mean(Y[,2]))/sqrt(var(Y[,2]))                      
Fr <- cca.results$V[,1:2]                                         
g <- cca.results$Gs[2,1:2]                                        

#points(cca.results$V[,1],cca.results$V[,2],cex=0.5,pch=19,col="red")
#textxy(cca.results$V[,1],cca.results$V[,2],rownames(Xn))

Cal.Data.Y2 <- calibrate(g,y,ticklabs,Fr,ticklab,lm=TRUE,tl=0.1,dp=TRUE,
                         reverse=TRUE,verb=TRUE,axislab="Y_2",
                         cex.axislab=0.75,showlabel=FALSE)

#cca.results<-lm.gls(Rxy[,5]~-1+Fr,W=solve(Rxx))



###################################################
### code chunk number 17: CalibrationGuide.Rnw:577-594
###################################################
data(linnerud)
X <- linnerud[,1:3]
Y <- linnerud[,4:6]
rda.results <- rda(X,Y)
plot(rda.results$Fs[,1],rda.results$Fs[,2],pch=16,asp=1,xlim=c(-2,2),ylim=c(-2,2),
     cex=0.5,xlab="1st principal axis",ylab="2nd principal axis")
arrows(0,0,2*rda.results$Gyp[,1],2*rda.results$Gyp[,2],length=0.1)
textxy(rda.results$Fs[,1],rda.results$Fs[,2],rownames(X),cex=0.75)
textxy(2*rda.results$Gyp[,1],2*rda.results$Gyp[,2],colnames(Y),cex=0.75)

y <- rda.results$Yh[,3]
g <- rda.results$Gyp[3,1:2]
Fr <- rda.results$Fs[,1:2]

ticklab <- c(seq(-0.6,-0.1,by=0.1),seq(0.1,0.6,by=0.1))
Calibrate.Yhat3 <- calibrate(g,y,ticklab,Fr,ticklab,lm=TRUE,dp=TRUE,tl=0.1,
                             axislab="Sauts",showlabel=FALSE)


###################################################
### code chunk number 18: CalibrationGuide.Rnw:597-621
###################################################
plot(rda.results$Gxs[,1],rda.results$Gxs[,2],pch=16,asp=1,xlim=c(-2,2),
     ylim=c(-2,2),cex=0.5,xlab="1st principal axis",
ylab="2nd principal axis")
arrows(0,0,rda.results$Gxs[,1],rda.results$Gxs[,2],length=0.1)
arrows(0,0,rda.results$Gyp[,1],rda.results$Gyp[,2],length=0.1)
textxy(rda.results$Gxs[,1],rda.results$Gxs[,2],colnames(X),cex=0.75)
textxy(rda.results$Gyp[,1],rda.results$Gyp[,2],colnames(Y),cex=0.75)

y <- rda.results$B[,3]
g <- rda.results$Gyp[3,1:2]
Fr <- rda.results$Gxs[,1:2]  

ticklab <- seq(-0.4,0.4,0.2)

W <-cor(X)

Calibrate.Y3 <- calibrate(g,y,ticklab,Fr,ticklab,lm=TRUE,dp=TRUE,tl=0.1,
                          weights=W,axislab="Sauts",showlabel=FALSE)
ticklab <- seq(-0.4,0.4,0.1)
Calibrate.Y3 <- calibrate(g,y,ticklab,Fr,ticklab,lm=FALSE,tl=0.05,verb=FALSE,
                          weights=W)
ticklab <- seq(-0.4,0.4,0.01)
Calibrate.Y3 <- calibrate(g,y,ticklab,Fr,ticklab,lm=FALSE,tl=0.025,verb=FALSE,
                          weights=W)


