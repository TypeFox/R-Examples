### R code from vignette source 'BsMD.Rnw'

###################################################
### code chunk number 1: BM86data
###################################################
options(width=80)
library(BsMD)
data(BM86.data,package="BsMD")
print(BM86.data)


###################################################
### code chunk number 2: BM86fitting
###################################################
advance.lm <- lm(y1 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    X10 + X11 + X12 + X13 + X14 + X15, data=BM86.data)
shrinkage.lm <- lm(y2 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    X10 + X11 + X12 + X13 + X14 + X15, data=BM86.data)
strength.lm <- lm(y3 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    X10 + X11 + X12 + X13 + X14 + X15, data=BM86.data)
yield.lm <- lm(y4 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
                    X10 + X11 + X12 + X13 + X14 + X15, data=BM86.data)
coef.tab <- data.frame(advance=coef(advance.lm),shrinkage=coef(shrinkage.lm),
                strength=coef(strength.lm),yield=coef(yield.lm))
print(round(coef.tab,2))


###################################################
### code chunk number 3: DanielPlots
###################################################
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.5,.5,0),oma=c(0,0,0,0),
    xpd=TRUE,pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
DanielPlot(advance.lm,cex.pch=0.8,main="a) Default Daniel Plot")
DanielPlot(advance.lm,cex.pch=0.8,main="b) Labelled Plot",pch=20,
    faclab=list(idx=c(2,4,8),lab=c(" 2"," 4"," 8")))


###################################################
### code chunk number 4: BsMD.Rnw:163-169
###################################################
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.5,.5,0),oma=c(0,0,0,0),
    xpd=TRUE,pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
DanielPlot(strength.lm,half=TRUE,cex.pch=0.8,main="a) Half-Normal Plot",
    faclab=list(idx=c(4,12,13),lab=c(" x4"," x12"," x13")))
DanielPlot(strength.lm,main="b) Normal Plot",
    faclab=list(idx=c(4,12,13),lab=c(" 4"," 12"," 13")))


###################################################
### code chunk number 5: BsMD.Rnw:207-216
###################################################
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(1.5,.5,0),oma=c(0,0,0,0),
    xpd=TRUE,pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
LenthPlot(shrinkage.lm)
title("a) Default Lenth Plot")
b <- coef(shrinkage.lm)[-1] #    Intercept removed
LenthPlot(shrinkage.lm,alpha=0.01,adj=0.2)
title(substitute("b) Lenth Plot (" *a* ")",list(a=quote(alpha==0.01))))
text(14,2*b[14],"P ",adj=1,cex=.7)  # Label x14 corresponding to factor P
text(15,2*b[15]," -M",adj=0,cex=.7) # Label x15 corresponding to factor -M


###################################################
### code chunk number 6: BsMD.Rnw:229-234
###################################################
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.5,.5,0),oma=c(0,0,0,0),
    pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
DanielPlot(yield.lm,cex.pch=0.6,main="a) Daniel Plot")
LenthPlot(yield.lm,alpha=0.05,xlab="factors",adj=.9,
main="b) Lenth Plot")


###################################################
### code chunk number 7: BsMD.Rnw:291-301
###################################################
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.5,.5,0),oma=c(0,0,0,0),
    pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
X <- as.matrix(BM86.data[,1:15])
y <- BM86.data[,16] # Using prior probability of p=0.20, and k=10 (gamma=2.49)
advance.BsProb <- BsProb(X=X,y=y,blk=0,mFac=15,mInt=1,p=0.20,g=2.49,ng=1,nMod=10)
print(advance.BsProb,X=FALSE,resp=FALSE,nMod=5)
plot(advance.BsProb,main="a) Bayes Plot")
DanielPlot(advance.lm,cex.pch=0.6,main="b) Daniel Plot",
        faclab=list(idx=c(2,4,8),lab=c(" x2"," x4"," x8")))
#title("Example I",outer=TRUE,line=-1,cex=.8)


###################################################
### code chunk number 8: BsMD.Rnw:324-337
###################################################
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.5,.5,0),oma=c(0,0,0,0),
    pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
X <- as.matrix(BM86.data[,1:15])
y <- BM86.data[,19]
# Using prior probability of p=0.20, and k=5,10,15
yield.BsProb <- BsProb(X=X,y=y,blk=0,mFac=15,mInt=1,p=0.20,g=c(1.22,3.74),ng=10,nMod=10)
summary(yield.BsProb)
plot(yield.BsProb,main="a) Bayes Plot")
#title(substitute("( " *g* " )",list(g=quote(1.2<=gamma<=3.7))),line=-1)
title(substitute("( " *g1* "" *g2* " )",list(g1=quote(1.2<=gamma),g2=quote(""<=3.7))),line=-1)
DanielPlot(yield.lm,cex.pch=0.6,main="b) Daniel Plot",
    faclab=list(idx=c(1,7,8,9,10,14),lab=paste(" ",c(1,7,8,9,10,14),sep="")))
#title("Example IV",outer=TRUE,line=-1,cex=.8)


###################################################
### code chunk number 9: BsMD.Rnw:367-385
###################################################
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.5,.5,0),oma=c(0,0,0,0),
    pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
data(BM93.e1.data,package="BsMD")
X <- as.matrix(BM93.e1.data[,2:6])
y <- BM93.e1.data[,7]
prob <- 0.25
gamma <- 1.6
# Using prior probability of p=0.20, and k=5,10,15
reactor5.BsProb <- BsProb(X=X,y=y,blk=0,mFac=5,mInt=3,p=prob,g=gamma,ng=1,nMod=10)
summary(reactor5.BsProb)
plot(reactor5.BsProb,main="a) Main Effects")

data(PB12Des,package="BsMD")
X <- as.matrix(PB12Des)
reactor11.BsProb <- BsProb(X=X,y=y,blk=0,mFac=11,mInt=3,p=prob,g=gamma,ng=1,nMod=10)
print(reactor11.BsProb,models=FALSE)
plot(reactor11.BsProb,main="b) All Contrasts")
#title("12-runs Plackett-Burman Design",outer=TRUE,line=-1,cex.main=0.9)


###################################################
### code chunk number 10: BsMD.Rnw:409-427
###################################################
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.5,.5,0),oma=c(0,0,1,0),
    pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
data(BM93.e2.data,package="BsMD")
X <- as.matrix(BM93.e2.data[,1:7])
y <- BM93.e2.data[,8]
prob <- 0.25
gamma <- c(1,2)
ng <- 20
# Using prior probability of p=0.20, and k=5,10,15
fatigueG.BsProb <- BsProb(X=X,y=y,blk=0,mFac=7,mInt=2,p=prob,g=gamma,ng=ng,nMod=10)
plot(fatigueG.BsProb$GAMMA,1/fatigueG.BsProb$prob[1,],type="o",
    xlab=expression(gamma),ylab=substitute("P{" *g* "|y}",list(g=quote(gamma))))
title(substitute("a) P{" *g* "|y}"%prop%"1/P{Null|y, " *g* "}",list(g=quote(gamma))),
    line=+.5,cex.main=0.8)
gamma <- 1.5
fatigue.BsProb <- BsProb(X=X,y=y,blk=0,mFac=7,mInt=2,p=prob,g=gamma,ng=1,nMod=10)
plot(fatigue.BsProb,main="b) Bayes Plot",code=FALSE)
title(substitute("( "*g*" )",list(g=quote(gamma==1.5))),line=-1)


###################################################
### code chunk number 11: BsMD.Rnw:453-470
###################################################
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(2,.5,0),oma=c(0,0,1,0),
    pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
data(BM93.e3.data,package="BsMD")
print(BM93.e3.data)
X <- as.matrix(BM93.e3.data[1:16,2:9])
y <- BM93.e3.data[1:16,10]
prob <- 0.25
gamma <- 2.0
# Using prior probability of p=0.25, and gamma=2.0
plot(BsProb(X=X,y=y,blk=0,mFac=8,mInt=3,p=prob,g=gamma,ng=1,nMod=10),
    code=FALSE,main="a) Fractional Factorial (FF)")

X <- as.matrix(BM93.e3.data[,c(2:9,1)])
y <- BM93.e3.data[,10]
plot(BsProb(X=X,y=y,blk=0,mFac=9,mInt=3,p=prob,g=gamma,ng=1,nMod=5),
    code=FALSE,main="b) FF with Extra Runs",prt=TRUE,)
mtext(side=1,"(Blocking factor blk)",cex=0.7,line=2.5)


###################################################
### code chunk number 12: MSBExample1
###################################################
par(mfrow=c(1,2),mar=c(3,4,1,1),mgp=c(2,.5,0),oma=c(0,0,1,0),
    pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
data(BM93.e3.data,package="BsMD")
X <- as.matrix(BM93.e3.data[1:16,c(1,2,4,6,9)])
y <- BM93.e3.data[1:16,10]
injection16.BsProb <- BsProb(X=X,y=y,blk=1,mFac=4,mInt=3,p=0.25,g=2,ng=1,nMod=5)

X <- as.matrix(BM93.e3.data[1:16,c(1,2,4,6,9)])
p <- injection16.BsProb$ptop
s2 <- injection16.BsProb$sigtop
nf <- injection16.BsProb$nftop
facs <- injection16.BsProb$jtop
nFDes <- 4
Xcand <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                    -1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,
                    -1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,
                    -1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,
                    -1,1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,-1),
                    nrow=16,dimnames=list(1:16,c("blk","A","C","E","H"))
                )
print(MD(X=X,y=y,nFac=4,nBlk=1,mInt=3,g=2,nMod=5,p=p,s2=s2,nf=nf,facs=facs,
            nFDes=4,Xcand=Xcand,mIter=20,nStart=25,top=5))


###################################################
### code chunk number 13: ReactorData
###################################################
data(Reactor.data,package="BsMD")
print(Reactor.data)
#print(cbind(run=1:16,Reactor.data[1:16,],run=17:32,Reactor.data[17:32,]))


###################################################
### code chunk number 14: BsMD.Rnw:572-595
###################################################
par(mfrow=c(1,2),mar=c(3,4,1,1),mgp=c(2,.5,0),oma=c(0,0,0,0),
    pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
fraction <- c(25,2,19,12,13,22,7,32)
cat("Fraction: ",fraction)
X <- as.matrix(cbind(blk=rep(-1,8),Reactor.data[fraction,1:5]))
y <- Reactor.data[fraction,6]
print(reactor8.BsProb <- BsProb(X=X,y=y,blk=1,mFac=5,mInt=3,
    p=0.25,g=0.40,ng=1,nMod=32),X=FALSE,resp=FALSE,factors=TRUE,models=FALSE)
plot(reactor8.BsProb,code=FALSE,main="a) Initial Design\n(8 runs)")
p <- reactor8.BsProb$ptop
s2 <- reactor8.BsProb$sigtop
nf <- reactor8.BsProb$nftop
facs <- reactor8.BsProb$jtop
nFDes <- 4
Xcand <- as.matrix(cbind(blk=rep(+1,32),Reactor.data[,1:5]))
print(MD(X=X,y=y,nFac=5,nBlk=1,mInt=3,g=0.40,nMod=32,p=p,s2=s2,nf=nf,facs=facs,
            nFDes=4,Xcand=Xcand,mIter=20,nStart=25,top=5),Xcand=FALSE,models=FALSE)
new.runs <- c(4,10,11,26)
cat("Follow-up:",new.runs)
X <- rbind(X,Xcand[new.runs,])
y <- c(y,Reactor.data[new.runs,6])
print(reactor12.BsProb <- BsProb(X=X,y=y,blk=1,mFac=5,mInt=3,p=0.25,g=1.20,ng=1,nMod=5))
plot(reactor12.BsProb,code=FALSE,main="b) Complete Design\n(12 runs)")


###################################################
### code chunk number 15: One-at-a-time
###################################################
data(Reactor.data,package="BsMD")
#cat("Initial Design:\n")
X <- as.matrix(cbind(blk=rep(-1,8),Reactor.data[fraction,1:5]))
y <- Reactor.data[fraction,6]
lst <- reactor8.BsProb <- BsProb(X=X,y=y,blk=1,mFac=5,mInt=3,p=0.25,g=0.40,ng=1,nMod=32)

#cat("Follow-Up: run 1\n")
p <- lst$ptop; s2 <- lst$sigtop; nf <- lst$nftop; facs <- lst$jtop
reactor8.MD <- MD(X=X,y=y,nFac=5,nBlk=1,mInt=3,g=0.40,nMod=32,p=p,s2=s2,nf=nf,facs=facs,
            nFDes=1,Xcand=Xcand,mIter=20,nStart=25,top=3)
new.run <- 10
X <- rbind(X,Xcand[new.run,]); rownames(X)[nrow(X)] <- new.run
y <- c(y,Reactor.data[new.run,6])
lst <- reactor9.BsProb <- BsProb(X=X,y=y,blk=1,mFac=5,mInt=3,p=0.25,g=0.7,ng=1,nMod=32)

#cat("Follow-Up: run 2\n")
p <- lst$ptop; s2 <- lst$sigtop; nf <- lst$nftop; facs <- lst$jtop
reactor9.MD <- MD(X=X,y=y,nFac=5,nBlk=1,mInt=3,g=0.7,nMod=32,p=p,s2=s2,nf=nf,facs=facs,
            nFDes=1,Xcand=Xcand,mIter=20,nStart=25,top=3)
new.run <- 4
X <- rbind(X,Xcand[new.run,]); rownames(X)[nrow(X)] <- new.run
y <- c(y,Reactor.data[new.run,6])
lst <- reactor10.BsProb <- BsProb(X=X,y=y,blk=1,mFac=5,mInt=3,p=0.25,g=1.0,ng=1,nMod=32)

#cat("Follow-Up: run 3\n")
p <- lst$ptop; s2 <- lst$sigtop; nf <- lst$nftop; facs <- lst$jtop
reactor10.MD <- MD(X=X,y=y,nFac=5,nBlk=1,mInt=3,g=1.0,nMod=32,p=p,s2=s2,nf=nf,facs=facs,
            nFDes=1,Xcand=Xcand,mIter=20,nStart=25,top=3)
new.run <- 11
X <- rbind(X,Xcand[new.run,]); rownames(X)[nrow(X)] <- new.run
y <- c(y,Reactor.data[new.run,6])
lst <- reactor11.BsProb <- BsProb(X=X,y=y,blk=1,mFac=5,mInt=3,p=0.25,g=1.3,ng=1,nMod=32)

#cat("Follow-Up: run 4\n")
p <- lst$ptop; s2 <- lst$sigtop; nf <- lst$nftop; facs <- lst$jtop
reactor10.MD <- MD(X=X,y=y,nFac=5,nBlk=1,mInt=3,g=1.3,nMod=32,p=p,s2=s2,nf=nf,facs=facs,
            nFDes=1,Xcand=Xcand,mIter=20,nStart=25,top=3)
new.run <- 15
X <- rbind(X,Xcand[new.run,]); rownames(X)[nrow(X)] <- new.run
y <- c(y,Reactor.data[new.run,6])
reactor12 <- BsProb(X=X,y=y,blk=1,mFac=5,mInt=3,p=0.25,g=1.30,ng=1,nMod=10)

print(reactor12,nMod=5,models=TRUE,plt=FALSE)


###################################################
### code chunk number 16: BsMD.Rnw:663-676
###################################################
par(mfrow=c(2,2),mar=c(3,4,1,1),mgp=c(2,.5,0),oma=c(1,0,1,0),
    pty="s",cex.axis=0.7,cex.lab=0.8,cex.main=0.9)
#plot(reactor8.BsProb,code=FALSE)
#mtext(side=1,"a) 8 runs",line=3,cex=0.7)
plot(reactor9.BsProb,code=FALSE)
mtext(side=1,"b) 9 runs",line=3,cex=0.7)
plot(reactor10.BsProb,code=FALSE)
mtext(side=1,"c) 10 runs",line=3,cex=0.7)
plot(reactor11.BsProb,code=FALSE)
mtext(side=1,"d) 11 runs",line=3,cex=0.7)
plot(reactor12.BsProb,code=FALSE)
mtext(side=1,"e) 12 runs",line=3,cex=0.7)
title("One-at-a-time Experiments",outer=TRUE)


