#Source: MOE-Jan21.R
#boostrap estimation of the MOE for the empirical standard devations
#
setwd("D:/DropBox/R/2015/artfima/simulations/Jan17")
FileSimulation<-"WS-Jan17-ll10000.RData"
attach(what=FileSimulation)
objects(2)


#simulation results#WS-Dec27-ll10000.RData
setwd("D:/DropBox/R/2015/artfima/simulations/Dec27")
FileSimulation<-"WS-Dec27-ll10000.RData"
attach(what=FileSimulation)
#check
objects(2)
#
#OUT components: rep(
#{n, lamda, d, lambdaHatWhittle, dHatWhittle, seWhittle}
# {n=c(200, 500, 1000,5000)})
#Note: se contains 2 elements, so there are 7 elements 
#
me <- matrix(0, ncol=9, nrow=NSIM*36)
i1 <- 1
for (i in 1:length(OUT)) {
  i2 <- i1+3
  me[i1:i2,] <- matrix(OUT[[i]], ncol=9, byrow=TRUE)
  i1 <- i2+1
}
me <- me[,-(8:9)]
dimnames(me)[[2]] <- names(me) <- c("n", "lambda", "d", "lambdaHat", "dHat", 
                                    "seLambdaHat", "seDHat")
medf <- as.data.frame(me) #large 360,000 rows
#
medf$n <- ordered(medf$n, levels=Ns, labels=paste0("n=",Ns))
medf$lambda <- ordered(medf$lambda, levels=LAMs, labels=paste0("lambda=",LAMs))
medf$d <- ordered(medf$d, levels=Ds, labels=paste0("d=",Ds))
#
#robust estimate of empirical sd's (use mad())
aLamHat <- with(medf, tapply(lambdaHat, list(n=n, lambda=lambda, d=d),mad))
adHat <- with(medf, tapply(dHat, list(n=n, lambda=lambda, d=d), mad))
vLamHat <- c(aLamHat)
vdHat <- c(adHat)
#
#Bootstrap
NBoot <- 1000 #155 sec
bmad <- function(x) {
  mad(sample(x, size=length(x), replace=TRUE))
}
startTime <- proc.time()[3]
set.seed(227771)
bdHat <- bLamHat <- matrix(numeric(0), nrow=NBoot, ncol=36)
for (iBoot in 1:NBoot) {
  bLamHat[iBoot,] <- c(with(medf,tapply(lambdaHat, list(n=n, lambda=lambda, d=d),
                                     bmad)))
  bdHat[iBoot,] <- c(with(medf, tapply(dHat, list(n=n, lambda=lambda, d=d),
                                    bmad)))
}
proc.time()[3]-startTime
#
boot.lam<-apply(bLamHat, 2, mad)
boot.d<-apply(bdHat, 2, mad)
layout(matrix(1:2, ncol=2))
boxplot(boot.lam, ylab=expression(paste("Bootstrap sd of empirical sd ",lambda)))
boxplot(boot.d, ylab="Bootstrap sd of empirical sd d")
layout(1)
summary(boot.lam)
summary(boot.d)
#
#table with empirical sd and their bootstrap 95% MOE
moe.lam <- 1.96*boot.lam
moe.d <- 1.96*boot.d
vLamHat <- c(aLamHat)
vdHat <- c(adHat)
vn <- ordered(rep(dimnames(aLamHat)[[1]],9), levels=dimnames(aLamHat)[[1]])
vLam <- ordered(rep(rep(dimnames(aLamHat)[[2]], rep(4, 3)),3), 
                levels=dimnames(aLamHat)[[2]])
vd <- ordered(rep(dimnames(aLamHat)[[3]], rep(12, 3)), 
              levels=dimnames(aLamHat)[[3]])
dfLam <- data.frame(n=vn, lam=vLam, d=vd, emsd.lamHat=vLamHat, emsd.d=vdHat, 
                    moe.lamHat=moe.lam, moe.d=moe.d)
#covert dfLam to nicer table
moe.tb <- dfLam
moe.tb$n <- as.integer(gsub(pattern="n=", replacement="", x=as.character(dfLam$n)))
moe.tb$lam <- as.numeric(gsub(pattern="lambda=", replacement="", x=as.character(dfLam$lam)))
moe.tb$d <- as.numeric(gsub(pattern="d=", replacement="", x=as.character(dfLam$d)))
moe.tb <- moe.tb[,c(1:4, 6, 5, 7)]
write.table(x=moe.tb, file="MOE.csv", row.names=FALSE)
#
#dotcharts with MOE
#lambda
MOE <- moe.tb$moe.lamHat
xyplot(factor(n)~emsd.lamHat | factor(d)*factor(lam), data=moe.tb,
       panel=function(x,y,subscripts=TRUE){
         panel.xyplot(x,y)
         panel.segments(x-MOE[subscripts], y, x+MOE[subscripts], y)
       },
       xlab=expression(sigma(hat(lambda))), ylab="n"
)
#d
MOE <- moe.tb$moe.d
xyplot(factor(n)~emsd.d | factor(d)*factor(lam), data=moe.tb,
       panel=function(x,y,subscripts=TRUE){
         panel.xyplot(x,y)
         panel.segments(x-MOE[subscripts], y, x+MOE[subscripts], y)
       },
       xlab=expression(sigma(hat(d))), ylab="n"
)
#
moe.tb












