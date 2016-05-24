## ----message=FALSE-------------------------------------------------------
library(lctools)
data(GR.Municipalities)
names(GR.Municipalities@data)

## ------------------------------------------------------------------------
Coords<-cbind(GR.Municipalities@data$X, GR.Municipalities@data$Y)
bw<-6
mI<-moransI(Coords,bw,GR.Municipalities@data$Income01)
moran.table<-matrix(data=NA,nrow=1,ncol=6)
colnames(moran.table) <- c("Moran's I", "Expected I", "Z resampling", "P-value resampling",
                     "Z randomization", "P-value randomization")
moran.table[1,1]<-mI$Morans.I
moran.table[1,2]<-mI$Expected.I
moran.table[1,3]<-mI$z.resampling
moran.table[1,4]<-mI$p.value.resampling
moran.table[1,5]<-mI$z.randomization
moran.table[1,6]<-mI$p.value.randomization

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(moran.table)

## ------------------------------------------------------------------------
bw<-c(3, 4, 6, 9, 12, 18, 24)

moran<-matrix(data=NA,nrow=7,ncol=7)
colnames(moran) <- c("ID", "k", "Moran's I", "Z resampling", "P-value resampling",
                     "Z randomization", "P-value randomization")
counter<-1

for(b in bw){
    moranI<-moransI(Coords,b,GR.Municipalities@data$Income01)
    moran[counter,1]<-counter
    moran[counter,2]<-b
    moran[counter,3]<-moranI$Morans.I
    moran[counter,4]<-moranI$z.resampling
    moran[counter,5]<-moranI$p.value.resampling
    moran[counter,6]<-moranI$z.randomization
    moran[counter,7]<-moranI$p.value.randomization
    counter<-counter+1
}

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(round(moran,4))

## ----fig.width = 5-------------------------------------------------------
plot(moran[,2], moran[,3], main="Global Moran's I", sub="", xlab="# of neighbours", 
     ylab="Moran's I")

## ------------------------------------------------------------------------
l.moran<-l.moransI(Coords,6,GR.Municipalities@data$Income01)

## ------------------------------------------------------------------------
xmin<-round(ifelse(abs(min(l.moran[,7])) > abs(min(l.moran[,8])), abs(min(l.moran[,7])), 
                   abs(min(l.moran[,8]))))
xmax<-round(ifelse(abs(max(l.moran[,7])) > abs(max(l.moran[,8])), abs(max(l.moran[,7])), 
                   abs(max(l.moran[,8]))))
xmax<-ifelse(xmin>xmax,xmin,xmax)+1
ymax<-xmax
xmin<- -xmax
ymin<- -ymax
reg1 <- lm(l.moran[,8]~l.moran[,7])

## ----fig.width = 5, fig.height = 5---------------------------------------
plot(l.moran[,7], l.moran[,8], main="Moran's I scatterplot", sub="", xlab="Income", 
     ylab="lagged Income", xlim=c(xmin, xmax), ylim=c(ymin, ymax))
abline(h=0)
abline(v=0)
abline(reg1, col="red")

