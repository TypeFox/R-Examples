#eSeriesB.R
#SeriesB dataset
#Exact MLE. Best Model.
#> TotalTimes
#   FD   FGN   PLA  NONE
# 5.19 17.64  5.71  1.68
#> sum(TotalTimes)
#[1]  30.22
#
require("FGN")
#z <- NileMin
z <- abs(diff(log(SeriesB)))
z <- z - mean(z)
n <- length(z)
P <- Q <- 3

TotalTimes <- numeric(4)
names(TotalTimes) <- c("FD", "FGN", "PLA", "NONE")
#FD
numMod <- (P+1)*(Q+1)
outMod <- vector("list", numMod)
ii <- 0
startTime <- proc.time()
for (p in 0:P)
    for (q in 0:Q) {
    ii <- ii+1
    k <- p+q+2
    order <- c(p,0,q)
    ans <- earfima(z, order=order, lmodel="FD")
    Le <- ans$LL
    bice <- -2*Le+k*log(n)
    out <- c(p,q,Le,bice)
    names(out) <- c("p","q","Le","bice")
    outMod[[ii]] <- out
    }
endTime <- proc.time()
totalTime <- endTime-startTime
TotalTimes[1] <- totalTime[1]
m<-matrix(unlist(outMod),byrow=TRUE,ncol=4)
dimnames(m)[[2]]<- c("p","q","Le","bice")
ind1 <- which.min(m[,"bice"])
mc<-rep(" ", 16)
mc[ind1]<-"*"
dimnames(m)[[1]]<-mc
mFD<-m
#

#FGN
numMod <- (P+1)*(Q+1)
outMod <- vector("list", numMod)
ii <- 0
startTime <- proc.time()
for (p in 0:P)
    for (q in 0:Q) {
    ii <- ii+1
    k <- p+q+2
    order <- c(p,0,q)
    ans <- earfima(z, order=order, lmodel="FGN")
    Le <- ans$LL
    bice <- -2*Le+k*log(n)
    out <- c(p,q,Le,bice)
    names(out) <- c("p","q","Le","bice")
    outMod[[ii]] <- out
    }
endTime <- proc.time()
totalTime <- endTime-startTime
TotalTimes[2] <- totalTime[1]
m<-matrix(unlist(outMod),byrow=TRUE,ncol=4)
dimnames(m)[[2]]<- c("p","q","Le","bice")
ind1 <- which.min(m[,"bice"])
mc<-rep(" ", 16)
mc[ind1]<-"*"
dimnames(m)[[1]]<-mc
mFGN<-m
#
#PLA
numMod <- (P+1)*(Q+1)
outMod <- vector("list", numMod)
ii <- 0
startTime <- proc.time()
for (p in 0:P)
    for (q in 0:Q) {
    ii <- ii+1
    k <- p+q+2
    order <- c(p,0,q)
    ans <- earfima(z, order=order, lmodel="PLA")
    Le <- ans$LL
    bice <- -2*Le+k*log(n)
    out <- c(p,q,Le,bice)
    names(out) <- c("p","q","Le","bice")
    outMod[[ii]] <- out
    }
endTime <- proc.time()
totalTime <- endTime-startTime
TotalTimes[3] <- totalTime[1]
m<-matrix(unlist(outMod),byrow=TRUE,ncol=4)
dimnames(m)[[2]]<- c("p","q","Le","bice")
ind1 <- which.min(m[,"bice"])
mc<-rep(" ", 16)
mc[ind1]<-"*"
dimnames(m)[[1]]<-mc
mPLA<-m
#
#
#NONE
numMod <- (P+1)*(Q+1)
outMod <- vector("list", numMod)
ii <- 0
startTime <- proc.time()
for (p in 0:P)
    for (q in 0:Q) {
    ii <- ii+1
    k <- p+q+2
    order <- c(p,0,q)
    ans <- earfima(z, order=order, lmodel="NONE")
    Le <- ans$LL
    bice <- -2*Le+k*log(n)
    out <- c(p,q,Le,bice)
    names(out) <- c("p","q","Le","bice")
    outMod[[ii]] <- out
    }
endTime <- proc.time()
totalTime <- endTime-startTime
TotalTimes[4] <- totalTime[1]
m<-matrix(unlist(outMod),byrow=TRUE,ncol=4)
dimnames(m)[[2]]<- c("p","q","Le","bice")
ind1 <- which.min(m[,"bice"])
mc<-rep(" ", 16)
mc[ind1]<-"*"
dimnames(m)[[1]]<-mc
mNONE<-m
#
LLs <- c(mFD["*",3],mFGN["*",3],mPLA["*",3],mNONE["*",3])
LLmax <- max(LLs)
RLs <- exp(LLs-LLmax)
names(RLs) <- c("FD","FGN","PLA","NONE")
#
bics <- c(mFD["*",4],mFGN["*",4],mPLA["*",4],mNONE["*",4])
bicmin <- min(bics)
RELs <- exp(-0.5*(bics-bicmin))
names(RELs) <- c("FD","FGN","PLA","NONE")
tb <- matrix(c(RLs,RELs)*100, byrow=TRUE, nrow=2)
dimnames(tb) <- list(c("RL","REL"), names(RELs))
tb
tbSeriesB <- tb
dump("tbSeriesB", file="d:/R/CRAN/FGN/vig/tbSeriesB.R")
#
TotalTimes
#
mFD
mFGN
mPLA
mNONE
