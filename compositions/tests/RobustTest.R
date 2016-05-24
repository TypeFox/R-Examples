
extensiveCheck <- FALSE
noRobustCheck  <- TRUE
if( !noRobustCheck ) {
###################################################
### chunk number 1: Bsp1
###################################################
#require(robustbase)
require(compositions)
#source("FunctionOutliers.R")
if( ! extensiveCheck ) {
simStore <-"storedsimulations.rda"
if( file.exists(simStore) )
  load(simStore,gsi.pStore)
}
#dades <- read.table("rutile.dat",header=TRUE)
#comp <- acomp(dades[,6:11]) 
#nom <- dades[,5] 
#location <- dades[,1] 
#gsize <-dades[,2] 
#losts <- dades[,13:18]
tmp<-set.seed(1400)
A <- matrix(c(0.1,0.2,0.3,0.1),nrow=2)
Mvar <- 0.1*ilrvar2clr(A%*%t(A))
Mcenter <- acomp(c(1,2,1))
typicalData <- rnorm.acomp(100,Mcenter,Mvar) # main population
colnames(typicalData)<-c("A","B","C")
data1 <- acomp(rnorm.acomp(100,Mcenter,Mvar))
data2 <- acomp(rbind(typicalData+rbinom(100,1,p=0.1)*rnorm(100)*acomp(c(4,1,1))))
data3 <- acomp(rbind(typicalData,acomp(c(0.5,1.5,2))))
colnames(data3)<-colnames(typicalData)
tmp<-set.seed(30)
# sapply(1:40,function(i) { set.seed(i);max(myMahalanobis.acomp(rnorm.acomp(100,Mcenter,Mvar)))})
rcauchy.acomp <- function (n, mean, var){
    D <- gsi.getD(mean)-1
    perturbe(ilrInv(matrix(rnorm(n*D)/rep(rnorm(n),D), ncol = D) %*% chol(clrvar2ilr(var))), mean)
}
data4 <- acomp(rcauchy.acomp(100,acomp(c(1,2,1)),Mvar/4))
colnames(data4)<-colnames(typicalData)
data5 <- acomp(rbind(unclass(typicalData)+outer(rbinom(100,1,p=0.1)*runif(100),c(0.1,1,2))))
data6 <- acomp(rbind(typicalData,rnorm.acomp(20,acomp(c(4,4,1)),Mvar)))
datas <- list(data1=data1,data2=data2,data3=data3,data4=data4,data5=data5,data6=data6)
tmp <-c()

if( !extensiveCheck )
  datas <- datas[c(2,5,6)]


###################################################
### chunk number 2: Fig1
###################################################

opar<-par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
tmp <- mapply(function(x,y) {
  plot(x,cex=0.5);
  title(y)
  meanRob<- mean(x,robust=TRUE)
  covRob <- var(x,robust=TRUE,giveCenter=TRUE)
  meanRob<- attr(covRob,"center")
  covCla <- cov(x)
  meanCla<- mean(x)
  ellipses(Mcenter,Mvar,col="green",r=3)
  ellipses(meanCla,covCla,col="red",r=3)
  ellipses(meanRob,covRob,col="blue",r=3)
},datas,names(datas))
tmp <- par(opar)


###################################################
### chunk number 3: Bsp eval=FALSE
###################################################
## outlierplot.acomp(comp,type="scatter",class.type="grade")


###################################################
### chunk number 4: Fig2aa
###################################################
opar<-par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
tmp<-mapply(function(x,y) {
outlierplot(x,type="scatter",class.type="grade");
  title(y)
},datas,names(datas))

#mapply(function(x,y) {
#  outlierplot.acomp(x,type="scatter",classifier=myClassifier3.acomp,class.type="best",legend=print(colcode));
#  title(y)
#},datas,names(datas))
tmp<-par(opar)


###################################################
### chunk number 5: Fig2
###################################################
opar<-par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
tmp<-mapply(function(x,y) {
  myCls2 <- OutlierClassifier1(x,alpha=0.05,type="all",corrected=TRUE)
  outlierplot(x,type="scatter",classifier=OutlierClassifier1,class.type="best",
  Legend=legend(1,1,levels(myCls),xjust=1,col=colcode,pch=pchcode),
  pch=as.numeric(myCls2));
  legend(0,1,legend=levels(myCls2),pch=1:length(levels(myCls2)))
  title(y)
},datas,names(datas))

#mapply(function(x,y) {
#  outlierplot.acomp(x,type="scatter",classifier=myClassifier3.acomp,class.type="best",legend=print(colcode));
#  title(y)
#},datas,names(datas))
tmp<-par(opar)


###################################################
### chunk number 6: Bsp eval=FALSE
###################################################
## outlierplot.acomp(comp,type="scatter",class.type="best")


###################################################
### chunk number 7: Bsp eval=FALSE
###################################################
## outlierplot.acomp(comp,type="ecdf")


###################################################
### chunk number 8: Bsp eval=FALSE
###################################################
## outlierplot.acomp(comp,type="nout")


###################################################
### chunk number 9: Fig3
###################################################
opar<-par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
for( i in 1:length(datas) ) 
  outlierplot(datas[[i]],type="ecdf",main=names(datas)[i])

if( FALSE ) {
for( i in 1:length(datas) ) {
  X <- datas[[i]]
  K <- sort(MahalanobisDist(X))
  Ks <- seq(min(K),max(K),length.out=100)
  cdf <- ecdf(K) 
  plot(cdf,main="")
  lines(Ks,pMahalanobis(Ks,N=nrow(X),d=gsi.getD(X)-1),col="red")
  KpQ <- pQuantileMahalanobis(Ks,N=nrow(X),d=gsi.getD(X)-1,p=0.05)
  lines(Ks,KpQ,col="green")
#  print(ifelse(cdf(K) < KpQ,1-KpQ/(1-cdf(K)),0))
  lines(Ks,ifelse(cdf(Ks) < KpQ,1-(1-KpQ)/(1-cdf(Ks)),0),col="gray") 
  title(main=names(datas)[i])
  invisible("")
}
} else if(FALSE)  {
for( i in 1:length(datas) ) {
  X <- datas[[i]]
  K <- sort(MahalanobisDist(X))
  Ks <- exp(seq(-5,5,length.out=100))
  cdf <- ecdf(K) 
  plot(cdf,main="",log="x",xlim=exp(c(-5,5)))
  lines(Ks,pMahalanobis(Ks,N=nrow(X),d=gsi.getD(X)-1),col="red",lwd=2)
  KpQ <- pQuantileMahalanobis(Ks,N=nrow(X),d=gsi.getD(X)-1,p=0.05)
  lines(Ks,KpQ,col="green",lwd=2)
#  print(ifelse(cdf(K) < KpQ,1-KpQ/(1-cdf(K)),0))
  lines(Ks,ifelse(cdf(Ks) < KpQ,1-(1-KpQ)/(1-cdf(Ks)),0),col="gray") 
  title(main=names(datas)[i])
  invisible("")
}
}


par(opar)


###################################################
### chunk number 10: Fig4aa
###################################################
opar<-par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1))  
for( i in 1:length(datas) )  {
  outlierplot(datas[[i]],type="nout",main=names(datas)[i])
  title(paste("data",i,sep=""))
}
if(FALSE) {
for( i in 1:length(datas) ) {
  X <- datas[[i]]
  K <- sort(myMahalanobis.acomp(X))
  n <- nrow(X)
  cdf <- ecdf(K) 
#  plot(cdf,main="")
  Kp  <- pMahalanobis(K,N=nrow(X),d=gsi.getD(X)-1)
#  lines(K,Kp,col="gray")
  KpQ  <- pQuantileMahalanobis(K,N=nrow(X),d=gsi.getD(X)-1,p=0.05)
  KpQ2 <- pQunatileMahalanobis(K,N=nrow(X),d=gsi.getD(X)-1,p=0.95)
  #lines(K,KpQ,col="gray",lty=2)
  #print(ifelse(cdf(K) < KpQ,1-KpQ/(1-cdf(K)),0))
  #print(KpQ2)
  plot(0,0,type="n",xlim=exp(c(-5,5)),ylim=n*range(c(Kp-cdf(K),Kp-KpQ,Kp-KpQ2)),log="x")
#  print(n*(Kp-KpQ2))
  lines(K,n*(Kp-cdf(K)),col="blue")
  lines(K,n*(Kp-KpQ),col="gray",lty=1)
#  lines(log(K),n*(Kp-KpQ2),col="green",lty=1)
  abline(h=0)
#  lines(K,ifelse(cdf(K) < KpQ,1-(1-KpQ)/(1-cdf(K)),0),col="red") 
  title(main=names(datas)[i])
}
}
par(opar)

if( extensiveCheck ) {

###################################################
### chunk number 11: Bsp2
###################################################
tmp<-set.seed(1400)
A <- matrix(c(0.1,0.2,0.3,0.1),nrow=2)
Mvar <- 0.1*ilrvar2clr(A%*%t(A))
Mcenter <- acomp(c(1,2,1))
typicalData <- rnorm.acomp(100,Mcenter,Mvar) # main population
colnames(typicalData)<-c("A","B","C")
data1a <- acomp(rnorm.acomp(100,Mcenter,Mvar))
colnames(data1a)<-colnames(typicalData)
data2a <- acomp(rbind(typicalData+rbinom(100,1,p=0.1)*rnorm(100)*0.4*acomp(c(4,1,1))))
data3a <- acomp(rbind(typicalData,acomp(c(0.1963,0.4618,0.3418))))
set.seed(30)
# sapply(1:40,function(i) { set.seed(i);max(myMahalanobis.acomp(rnorm.acomp(100,Mcenter,Mvar)))})
rmixgauss.acomp <- function (n, mean, var){
    D <- gsi.getD(mean)-1
    perturbe(ilrInv(matrix(rnorm(n*D)*rep(rexp(n),D), ncol = D) %*% chol(clrvar2ilr(var))), mean)
}
data4a <- acomp(rmixgauss.acomp(100,acomp(c(1,2,1)),Mvar))
colnames(data4a)<-colnames(typicalData)
data5a <- acomp(rbind(unclass(typicalData)+outer(0.1*rbinom(100,1,p=0.1)*runif(100),c(0.1,1,2))))
data6a <- acomp(rbind(typicalData,rnorm.acomp(20,acomp(c(0.28,0.49,0.21)),Mvar)))
datasA <- list(data1a=data1a,data2a=data2a,data3a=data3a,data4a=data4a,data5a=data5a,data6a=data6a)
tmp<-c()

}
###################################################
### chunk number 12: FigCluDir
###################################################
par(mfrow=c(2,3))
takeIf <- function(c,x,y) {
    y <- ifelse(c,y,y)
    y[c]<-x
    y
  }
if( extensiveCheck ) {

mapply(function(x,y) {
  center <- mean(x,robust=TRUE)
  myCls <- OutlierClassifier1(x,type="outlier")
  out <- myCls!="ok"
  if( sum(out+0)>2 ) {
    d <- asin(dist(normalize(acomp(x[out,])-center))/2)*2*180/pi
    hc <- hclust(d,"complete")
    trx <- cutree(hc,h=15)
    cols <- takeIf(out,as.numeric(trx)+1,1)
    plot(acomp(x),col=cols,pch=cols,main=y)
    title(paste(y,"Aitchison geometry",sep="\n",col="\n"))
  } else {
    plot.new()
  }
},datas[c(2,5,6)],names(datas)[c(2,5,6)])
mapply(function(x,y) {
  center <- mean(rcomp(x))
  myCls <- OutlierClassifier1.acomp(x,type="outlier")
  out <- myCls!="ok"
  if( sum(out+0)>2 ) {
    d <- asin(dist(normalize(ipt(rcomp(x[out,]))-ipt(center)))/2)*2*180/pi
    hc <- hclust(d,"complete")
    trx <- cutree(hc,h=15)
    cols <- takeIf(out,as.numeric(trx)+1,1)
    plot(rcomp(x),col=cols,pch=cols,main=y)
    title(paste(y,"real geometry",sep="\n",col="\n"))
  } else {
    plot.new()
  }
},datas[c(2,5,6)],names(datas)[c(2,5,6)])

# The following is much to time expensive for checking
#
}
if( extensiveCheck ) {

###################################################
### chunk number 13: ReadRutile
###################################################
#dades <- read.table("rutile.dat",header=TRUE)
data(Hydrochem)
comp <- acomp(Hydrochem[,6:12])
location <- Hydrochem[,"Location"]
gsize    <- Hydrochem[,"River"]
row.names(comp)<- paste("A",1:nrow(comp))

###################################################
### chunk number 14: FigRutType1
###################################################
 outlierplot(comp,type="scatter",class.type="grade",main="" )
#,
#  Legend=PSlegend("FigRutileType1Leg.eps",
#    legend=levels(myCls),col=colcode,pch=pchcode,horizontal=TRUE,ncol=length(levels(myCls)))


###################################################
### chunk number 15: FigRutType3
###################################################
 par(mfrow=c(1,2),mar=c(4,4,1,1))
 outlierplot(comp,type="ecdf")
 outlierplot(comp,type="nout")


###################################################
### chunk number 16: FigRutType2
###################################################
 outlierplot.acomp(comp,type="scatter",class.type="best",center=TRUE)


###################################################
### chunk number 17: FigRutClu1
###################################################
  x<-comp
  center <- mean(rcomp(x),robust=TRUE)
  myCls <- OutlierClassifier1(x,type="outlier")
  out <- myCls!="ok"
  d <- asin(dist(normalize(ipt(rcomp(x[out,]))-ipt(center)))/2)*2*180/pi
  hc <- hclust(d,"complete")
  par(mar=c(1,4,0,0))
  plot(hc,main="",xlab="")


###################################################
### chunk number 18: FigRutClu1a
###################################################
  trx1 <- cutree(hc,h=50)
  cols <- colorsForOutliers1(factor(trx1))
  plot(acomp(x),col=takeIf(out,cols[as.numeric(trx1)],"gray"),
   pch=takeIf(out,as.numeric(trx1)+1,"."),center=TRUE)


###################################################
### chunk number 19: FigRutClu2
###################################################
  x<-comp
  center <- mean(rcomp(x),robust=TRUE)
  myCls <- OutlierClassifier1(x,type="outlier")
  out <- myCls!="ok"
  d <- dist(acomp(x[out,]))
  hc <- hclust(d,"complete")
  par(mar=c(1,4,0,0))
  plot(hc,main="",xlab="")


###################################################
### chunk number 20: FigRutClu2a
###################################################
  trx2 <- cutree(hc,k=5)
  cols <- colorsForOutliers1(factor(trx2))
  plot(acomp(x),col=takeIf(out,cols[as.numeric(trx2)],"gray"),
       pch=takeIf(out,as.numeric(trx2)+1,"."),center=TRUE)


###################################################
### chunk number 21: FigRutileGlobal2
###################################################
  cl2<- OutlierClassifier1(comp,alpha=0.05,type="best")
par(mfrow=c(1,2))
 trxext1 = as.integer(cl2)*0+1
   names(trxext1)=row.names(comp)
   trxext1[names(trx1)]=trx1+1
  cols <- colorsForOutliers1(cl2)
  barplot(table(cl2,trxext1),col=cols,main="")
  legend(x=3,y=120, legend=levels(cl2),fill=cols)
 trxext2 = as.integer(cl2)*0+1
   names(trxext2)=row.names(comp)
   trxext2[names(trx2)]=trx2+1
  cols <- colorsForOutliers1(cl2)
  barplot(table(cl2,trxext2),col=cols,main="")
#  legend(x=0,y=1,
#      legend=levels(cl2),col=colcode,pch=pchcode,horizontal=TRUE,
#      ncol=length(levels(myCls)))


###################################################
### chunk number 22: FigRutileGlobal2a
###################################################
 outlierplot(comp,type="biplot",class.type="best",main="",pch=cl$types,
  col=colorsForOutliers1) 


###################################################
### chunk number 23: FigRutileRobBiplot
###################################################
      colcode <- colorsForOutliers1(cl2)
      pchcode <- pchForOutliers1(cl2)
         pchcode[1]='\020'
      pch <- pchcode[as.integer(cl2)]
      xcol <- colcode[as.integer(cl2)]

      coloredBiplot(xrf=svd(var(comp,robust=TRUE)),xcol=xcol,pc.biplot=FALSE,pch=pch)
var(comp,robust=TRUE,giveCenter=TRUE)

###################################################
### chunk number 24: FigRutileSubcomp
###################################################
par(mfrow=c(2,2),oma=c(0,0,0,0))

xcol = ifelse(xcol=="gray40","white",xcol)
 plot(acomp(comp[,1:3]),bg=xcol,pch=c(1:3,21,5,24)[as.integer(location)],
     center=TRUE,cex=1.25)
 legend(x=0,y=1,legend=levels(location),pch=c(1:3,21,5,24))

 plot(acomp(comp[,1:3]),bg=xcol,pch=c(21,24,5)[as.integer(gsize)],center=TRUE,cex=1.25)
 legend(x=0,y=1,legend=levels(gsize),pch=c(21,24,5))

 options(robust=TRUE)
 plot(acomp(comp[,1:3]),bg=xcol,pch=c(1:3,21,5,24)[as.integer(location)],
     center=TRUE,cex=1.25)
 legend(x=0,y=1,legend=levels(location),pch=c(1:3,21,5,24))

 plot(acomp(comp[,1:3]),bg=xcol,pch=c(21,24,5)[as.integer(gsize)],center=TRUE,cex=1.25)
 legend(x=0,y=1,legend=levels(gsize),pch=c(21,24,5))
 options(robust=FALSE)

############# End of rutile figures
}

###################################################
### chunk number 25: Bsp eval=FALSE
###################################################
## outlierplot.acomp(comp, type, class.type, ...)


###################################################
### chunk number 26: FigXX2
###################################################
opar<-par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1)) 
tmp<-mapply(function(x,y) { 
  cl <- ClusterFinder1(x,sigma=0.4,radius=1) 
  plot(x,col=as.numeric(cl$types),pch=as.numeric(cl$types))
  legend(1,1,legend=levels(cl$types),xjust=1,col=1:length(levels(cl$types)),pch=1:length(levels(cl$types)))
  title(y)
  },datas,names(datas)) 
tmp<-par(opar)


###################################################
### chunk number 27: FigXX4
###################################################
opar<-par(mfrow=c(2,3),pch=19,mar=c(3,2,2,1)) 
tmp<-mapply(function(x,y) { 
  cl <- ClusterFinder1(x,sigma=0.4,radius=1) 
  cl2<- OutlierClassifier1(x,alpha=0.05,type="best")
  cols <- 1:length(levels(cl2))
  barplot(table(cl2,cl$types),col=cols,main=paste("Groups in",y),ylim=c(0,100)) 
  if( length(levels(cl$types))>1 )
    legend(par("usr")[2],100,xjust=1,levels(cl2),fill=cols)
  },datas,names(datas)) 
tmp<-par(opar)


###################################################
### chunk number 28: Bsp1 eval=FALSE
###################################################
## require(rrcov)
## require(compositions)
## source("FunctionOutliers.R")
## # covariance matrix
## A <- matrix(c(0.1,0.2,0.3,0.1),nrow=2)
## Mvar <- 0.1*ilrvar2clr(A%*%t(A))
## # mean composition
## Mcenter <- acomp(c(1,2,1))
## # basic simulation
## typicalData <- rnorm.acomp(100,Mcenter,Mvar) # main population
## colnames(typicalData)<-c("A","B","C")
## 
## # data1: without outliers
## data1 <- acomp(rnorm.acomp(100,Mcenter,Mvar))
##  colnames(data1)<-colnames(typicalData)
## 
## # data2: with 10% of data with an error in component A
## data2 <- acomp( rbind( typicalData + 
##        rbinom(100,1,p=0.1)*rnorm(100)*acomp(c(4,1,1)) )  )
## 
## # data3: with an erratic outlier
## data3 <- acomp(rbind(typicalData,acomp(c(0.5,1.5,2))))
## 
## # data4: with heavy tails (Cauchy type)
## tmp<-set.seed(30)
## rcauchy.acomp <- function (n, mean, var){
##  D <- gsi.getD(mean)-1
##  perturbe( 
##    ilrInv(
##     matrix(rnorm(n*D)/rep(rnorm(n),D), ncol = D) %*% chol(clrvar2ilr(var))
##    ), mean)
## }
## data4 <- acomp(rcauchy.acomp(100,acomp(c(1,2,1)),Mvar/4))
## colnames(data4)<-colnames(typicalData)
## 
## # data5: with an additive pollution through [A,B,C]=c(0.1,1,2)
## data5 <- acomp( rbind( unclass(typicalData)+
##        outer(rbinom(100,1,p=0.1)*runif(100),c(0.1,1,2)) ) )
## 
## # data6: with two subgroups: typicalData,
## #     and a newly generated one, around [A,B,C]=c(4,4,1)
## data6 <- acomp( rbind(
##           typicalData,
##           rnorm.acomp(20,acomp(c(4,4,1)),Mvar) ) )

if( !extensiveCheck )
  try(save(file=simStore,list=objects(gsi.pStore),envir=gsi.pStore))
}
