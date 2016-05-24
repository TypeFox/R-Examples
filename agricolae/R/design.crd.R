`design.crd` <-
function(trt,r,serie=2,seed=0,kinds="Super-Duper",randomization=TRUE)
{
number<-0
if(serie>0) number<-10^serie
junto<-data.frame(trt,r)
junto<-junto[order(junto[,1]),]
TR<-as.character(junto[,1])
r<-as.numeric(junto[,2])
y <- rep(TR[1], r[1])
tr <- length(TR)
if (seed == 0) {
genera<-runif(1)
seed <-.Random.seed[3]
}
set.seed(seed,kinds)
parameters<-list(design="crd",trt=trt,r=r,serie=serie,seed=seed,kinds=kinds,randomization)
for (i in 2:tr) y <- c(y, rep(TR[i], r[i]))
trat<-y
if(randomization)trat <- sample(y, length(y), replace = FALSE)
	plots <- number+1:length(trat)
	dca<-data.frame(plots, trat)
	dca[,1]<-as.numeric(dca[,1])
	xx<-dca[order(dca[,2],dca[,1]),]
	r1<-seq(1,r[1])
for (i in 2:length(r)) {
	r1<-c(r1,seq(1,r[i]))
}
yy<-data.frame(plots=xx[,1],r=r1,xx[,2])
book<-yy[order(yy[,1]),]
rownames(book)<-rownames(yy)
names(book)[3]<-c(paste(deparse(substitute(trt))))
outdesign<-list(parameters=parameters,book=book)
return(outdesign)
}

