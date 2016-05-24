#Source: SimAnalyze-Jan19.R
require("lattice")
require("plyr") #need for maply()
informationMatrixARTFIMA <- function(lambda, d){
  v11 <- 4*pi*exp(-2*lambda)*(1+ 0.5*exp(-2*lambda)) #eqn (43)
  v22 <- (d^2)*exp(-2*lambda)/(1-exp(-2*lambda)) #eqn (48)
  v12 <- d*log(1-exp(-2*lambda)) #eqn (54)
  matrix(c(v11,v12,v12,v22), ncol=2, 
         dimnames=list(c("lambda","d"),c("lambda","d")))
}
seARTFIMA <- function(lambda, d, n){
  v <- solve(informationMatrixARTFIMA(lambda, d)*n)
  sqrt(diag(v))
}
#
#analysis of simulation results
setwd("D:/DropBox/R/2015/artfima/simulations/Jan17")
FileSimulation<-"WS-Jan17-ll10000.RData"
attach(what=FileSimulation)
#check
objects(2)
NSIM
LAMs
Ds
Ns
length(OUT)
OUT[[1]][1:9]
OUT[[1]]
length(OUT[[1]]/9)
#
#OUT components: rep(
#{n, lamda, d, lambdaHatWhittle, dHatWhittle, seWhittle}
# {n=c(200, 500, 1000,5000)})
#Note: se contains 2 elements, so there are 7 elements 
#
#check converge/boundary issues
me <- matrix(0, ncol=9, nrow=NSIM*36)
i1 <- 1
for (i in 1:length(OUT)) {
  i2 <- i1+3
  me[i1:i2,] <- matrix(OUT[[i]], ncol=9, byrow=TRUE)
  i1 <- i2+1
}
dimnames(me)[[2]] <- names(me) <- c("n", "lambda", "d", "lambdaHat", "dHat", 
                      "seLambdaHat", "seDHat", "convQ", "boundQ")
#
#check for convergence issues
table(me[,"convQ"])
table(me[,"boundQ"])
indBound <- me[,"boundQ"]>0
meB <- me[indBound,]
meB <- data.frame(meB)
with(meB, tapply(boundQ, list(n,lambda, d),  sum))
#
me <- me[,1:7]
medf <- as.data.frame(me) #large 360,000 rows
medf$n <- ordered(medf$n, levels=Ns, labels=paste0("n=",Ns))
medf$lambda <- ordered(medf$lambda, levels=LAMs, labels=paste0("lambda=",LAMs))
medf$d <- ordered(medf$d, levels=Ds, labels=paste0("d=",Ds))
#
#distribution of lambda-hat is highly skewed with outliers when n=200 improves
# as n increases so when n=5000 much better.
bwplot(n~lambdaHat|d*lambda, data=medf,layout=c(3,3), xlab=expression(hat(lambda)),
       main=expression("Box plots based on 10,000 replications"),
       panel=function(x,y){
         panel.bwplot(x,y)
         panel.grid(v=-1)
       }, scales=list(x=list(relation="free")))
#drilling down to n=5000
bwplot(n~lambdaHat|d*lambda, data=medf,layout=c(3,3), xlab=expression(hat(lambda)),
       main=expression("Box plots based on 10,000 replications"),
       panel=function(x,y){
         panel.bwplot(x,y)
         panel.grid(v=-1)
       }, scales=list(x=list(relation="free")),
       subset = n=="n=5000")
#
#distribution of d-hat
bwplot(n~dHat|lambda*d, data=medf,layout=c(3,3), xlab="d",
       main=expression("Box plots based on 10,000 replications"),
       panel=function(x,y){
         panel.bwplot(x,y)
         panel.grid(v=-1)
       }, scales=list(x=list(relation="free")))
bwplot(n~dHat|lambda*d, data=medf,layout=c(3,3), xlab="d",
       main=expression("Box plots based on 10,000 replications"),
       subset = n=="n=5000",
       panel=function(x,y){
         panel.bwplot(x,y)
         panel.grid(v=-1)
       }, scales=list(x=list(relation="free")))
#robust estimate of empirical sd's (use mad())
aLamHat <- with(medf, tapply(lambdaHat, list(n=n, lambda=lambda, d=d),mad)) 
adHat <- with(medf, tapply(dHat, list(n=n, lambda=lambda, d=d), mad)) 
#robust estimate of observed sd's (use median())
obseLamHat <- with(medf, tapply(seLambdaHat, list(n=n, lambda=lambda, d=d),median))  
obsedHat <- with(medf, tapply(seDHat, list(n=n, lambda=lambda, d=d), median)) 
#vectorize result and make data.frame, dfLam
vLamHat <- c(aLamHat)
vdHat <- c(adHat)
vobseLamHat <- c(obseLamHat)
vobsedHat <- c(obsedHat)
vn <- ordered(rep(dimnames(aLamHat)[[1]],9), levels=dimnames(aLamHat)[[1]])
vLam <- ordered(rep(rep(dimnames(aLamHat)[[2]], rep(4, 3)),3), 
                levels=dimnames(aLamHat)[[2]])
vd <- ordered(rep(dimnames(aLamHat)[[3]], rep(12, 3)), 
              levels=dimnames(aLamHat)[[3]])
dfLam <- data.frame(emsd.lamHat=vLamHat, emsd.d=vdHat, n=vn, lam=vLam, d=vd,
                    obse.lamHat=vobseLamHat, obse.dHat=vobsedHat)
#compute se with seARTFIMA
ind <- c(t(outer(seq(0, 320000, by=4*10^4), 1:4, "+")))
out <- maply(me[ind,c(2,3,1)], seARTFIMA)
seLam <- c(aperm(out[,,,1], c(3,1,2)))
sed <- c(aperm(out[,,,2], c(3,1,2)))
dfLam2 <- cbind(dfLam, thse.lam=seLam, thse.d=sed)
#barchart and dotchart comparison
barchart(d~emsd.lamHat|n*lam, data=dfLam2)
dotplot(d~emsd.lamHat|n*lam, data=dfLam2)
#dotplot with groups
#first make new data frame. combine emsd & thse
outLam <- make.groups(dfLam2$emsd.lamHat, dfLam2$thse.lam)
levels(outLam$which)<-c("emse","asyse")
outD <- make.groups(dfLam2$emsd.d, dfLam2$thse.d)
levels(outD$which)<-c("emse","asyse")
names(outLam)[1] <- "lambda"
names(outD)[1] <- "d"
out3 <- cbind(outLam, outD) #df with c(em,th) for lam & d
dimnames(out3)[[1]] <- 1:nrow(out3)
out3 <- out3[,-4] #out3 is 72-by-3 
#add ob both lam & d, creating 108-by-3 out4
out4lambda <- c(out3$lambda, dfLam2$obse.lamHat)
out4d <- c(out3$d, dfLam2$obse.dHat)
out4which <- ordered(c(as.character(out3$which), rep("obse", 36)))
out4 <- data.frame(lambda=out4lambda, d=out4d, which=out4which)
fn <- rep(dfLam2$n, 3)
flam <- rep(dfLam2$lam, 3)
fd <- rep(dfLam2$d, 3)
dfLam3 <- cbind(out4, fn=fn, flam=flam, fd=fd)
test <- ordered(dfLam3$which, levels=c("emse", "obse", "asyse"))
dfLam3$which <- test
#
#now dotplots - comparing observed and expected se methods
#lambda
dotplot(fn~lambda | fd*flam, groups=which, data=dfLam3, pch=c(8,15,19), 
        cex=2.5, col=c("brown",rgb(0,1,0,0.6),rgb(1,0,0,0.6)),
        xlab=expression(sigma(hat(lambda))),
        key=list(points=list(pch=c(8,15,19), col=c("brown",rgb(0,1,0,1),rgb(1,0,0,1))),
                 text=list(lab=c("empirical","observed","expected")),
                 columns=3,
                 title="Comparison of Empirical Sd, Observed & Expected SE"),
        )

#d
dotplot(fn~d | fd*flam, groups=which, data=dfLam3, pch=c(8,15,19), 
        cex=2.5, col=c("brown",rgb(0,1,0,0.6),rgb(1,0,0,0.6)), 
        xlab=expression(sigma(hat(d))),
        key=list(points=list(pch=c(8,15,19), col=c("brown",rgb(0,1,0,1),rgb(1,0,0,1))),
                 text=list(lab=c("empirical","observed","asymptotic")),
                 columns=3,
                 title="Comparison of Empirical Sd, Observed & Asymptotic SE")
)
#
#dotcharts with true asymptotic rather than estimated
n <- as.integer(gsub(pattern="n=", replacement="", x=as.character(dfLam$n)))
lam <- as.numeric(gsub(pattern="lambda=", replacement="", x=as.character(dfLam$lam)))
d <- as.numeric(gsub(pattern="d=", replacement="", x=as.character(dfLam$d)))
asySDs <- maply(cbind(lambda=lam[1:36], d=d[1:36], n=n[1:36]), seARTFIMA, .expand=FALSE)
dfLam4 <- dfLam3 
dfLam4[37:72,1:2] <- asySDs
dfLam4$alambda <- sqrt(n)*dfLam4$lambda
dfLam4$ad <- sqrt(n)*dfLam4$d
#un-adjusted se
dotplot(flam~lambda | fd*fn, groups=which, data=dfLam4, pch=c(8,15,19), 
        cex=2.5, col=c("brown",rgb(0,1,0,0.6),rgb(1,0,0,0.6)),
        xlab=expression(sigma(hat(lambda))),
        key=list(points=list(pch=c(8,15,19), col=c("brown",rgb(0,1,0,1),rgb(1,0,0,1))),
                 text=list(lab=c("empirical","observed","asymptotic")),
                 columns=3,
                 title="Comparison of Empirical Sd, Observed & Asymptotic SE"),
        scales=list(x=list(relation="free", tick.number=3))
)
#Revised plots##################################################################
#se's have been standardized as well
#
#lambda: empirical, observed, asymptotic
#1: raw. Better reflects the convergence with n perhaps.
dotplot(fn~lambda | fd*flam, groups=which, data=dfLam4,
        xlab=expression(sigma(hat(lambda))), pch=c(8,15,19),cex=c(2.5,2.5,1.2),
        col=c("brown",rgb(0,1,0,0.6),rgb(0,0,0,1)),
        key=list(points=list(pch=c(8,15,19), col=c("brown",rgb(0,1,0,1),rgb(0,0,0,1))),
                 text=list(lab=c("empirical","observed","asymptotic")),columns=3,
                 title="Empirical Sd, Observed & Asymptotic SE")
)

#2: standardized. In some cases, convergence is apparent but less so in others
dotplot(fn~alambda | fd*flam, groups=which, data=dfLam4,
        xlab=expression(sigma(hat(lambda))), pch=c(8,15,19),cex=c(2.5,2.5,1.2),
        col=c("brown",rgb(0,1,0,0.6),rgb(0,0,0,1)),
        key=list(points=list(pch=c(8,15,19), col=c("brown",rgb(0,1,0,1),rgb(0,0,0,1))),
                 text=list(lab=c("empirical","observed","asymptotic")),columns=3,
                 title="Standardized Empirical Sd, Observed & Asymptotic SE")
)

#d: empirical, observed, asymptotic
#1: raw. Better reflects the convergence with n perhaps.
dotplot(fn~d | fd*flam, groups=which, data=dfLam4,
        xlab=expression(sigma(hat(d))), pch=c(8,15,19),cex=c(2.5,2.5,1.2),
        col=c("brown",rgb(0,1,0,0.6),rgb(0,0,0,1)),
        key=list(points=list(pch=c(8,15,19), col=c("brown",rgb(0,1,0,1),rgb(0,0,0,1))),
                 text=list(lab=c("empirical","observed","asymptotic")),columns=3,
                 title="Empirical Sd, Observed & Asymptotic SE")
)

#2: standardized. In some cases, convergence is apparent but less so in others
dotplot(fn~ad | fd*flam, groups=which, data=dfLam4,
        xlab=expression(sigma(hat(d))), pch=c(8,15,19),cex=c(2.5,2.5,1.2),
        col=c("brown",rgb(0,1,0,0.6),rgb(0,0,0,1)),
        key=list(points=list(pch=c(8,15,19), col=c("brown",rgb(0,1,0,1),rgb(0,0,0,1))),
                 text=list(lab=c("empirical","observed","asymptotic")),columns=3,
                 title="Standardized Empirical Sd, Observed & Asymptotic SE")
)

#generate se.csv file, uncomment code below
#tb <- dfLam2[,c(1,2,6,7,8,9)]
tb <- dfLam2
tb$n <- as.integer(gsub(pattern="n=", replacement="", x=as.character(dfLam$n)))
tb$lam <- as.numeric(gsub(pattern="lambda=", replacement="", x=as.character(dfLam$lam)))
tb$d <- as.numeric(gsub(pattern="d=", replacement="", x=as.character(dfLam$d)))
tb <- tb[,c(3,4,5, 1, 2, 6, 7, 8, 9)]
names(tb) <- gsub(pattern="th", replacement="ex", names(tb))
names(tb) <- gsub(pattern="Hat", replacement="", names(tb))
write.table(x=tb, file="se.csv", row.names=FALSE)
tbasy <- cbind(tb[1:36,1:3], asySDs)
write.table(x=tb, file="asyse.csv", row.names=FALSE)
#
#raw R values
print(tb, digits=4)
print(tbasy, digits=4)

