stkde.sig<-function(xlong,ylat,ztime,xgrids,ygrids,breaks=0.05,sim=100,alpha=0.05,nrowspar=1,...) {
  #=======================================================
  #
  #  TITLE:     Spatio-Temporal Kernel Density Estimation with Significant P-Value contours
  #  FUNCTION:  stkde.sig()
  #  AUTHOR:    Zhijie Zhang
  #  DATE:      20 JANUARY 2010
  #  CALLS:     np,sp,graphics
  #  NEEDS:
  #  NOTES:
  #  CHANGE HISTORY:
  #=======================================================
# SET DEPENDENCIES
require(np)
require(sp)
#create necessary datasets for intermediate manipulations
tlength <- length(sort(unique(ztime)))
x.seq0 <- seq(floor(min(xlong)),ceiling(max(xlong)),length=xgrids)
y.seq0 <- seq(floor(min(ylat)),ceiling(max(ylat)),length=ygrids)

f<-list()
Myregion <- SpatialPoints(cbind(xlong,ylat))
samsize<-length(xlong)
prank<-array(NA,c(xgrids,ygrids,tlength))

#my data results 
#f[[1]] <- stkde.base(xlong,ylat,ztime,xgrids=20,ygrids=20,breaks=0.05)
f[[1]] <- stkde.base(xlong,ylat,ztime,xgrids,ygrids,breaks)
#Simulated results
for (nsim in 2:sim)  {
samdata<-spsample(Myregion,samsize,type="random",bb=bbox(Myregion))
samdata<-data.frame(samdata)
#divide "samdata" into "tlength" parts
#split(dataset,cut(1:nrow(dataset),5))
samdata$ztime <- cut(samdata[,1], breaks=tlength, labels =c(1:tlength),ordered = TRUE)
#colnames(samdata) <- c("xlongsim","ylatsim","ztimesim")
attach(samdata)
f[[nsim]] <- stkde.base(samdata[,1],samdata[,2],samdata[,3],xgrids,ygrids,breaks)
detach(samdata)
}
#Now obtain the ranks for the my dataset,f[[1]]
#in f[[1]], and all the simulated f[[2]]--f[[nsim]] results
for (i in 1:xgrids) {
 for (j in 1:ygrids) {
  for (k in 1:tlength) {
     prank[i, j, k]<-rank(unlist(lapply(f, '[', i, j, k)))[1]
    }
  }
}
#Caculating the pvalue with prank
pvalue<-(prank+1)/sim
#transform the "pvalue" array into dataframe "pvalues"
#pvalues<-sapply(as.data.frame.table(pvalue), as.numeric)
#sort.aav<-aav[order(aav$aaa,aav$aaa2),]
#Now we get the arrays,"f[[1]]"  and "pvalue", for our dataset, kde and P values.
#Now we can map the "KDE risk map" and "p value surface"
if (tlength%%2==0)  {
 #mapping the results,even number
 brks <- quantile(f[[1]], seq(0,1,breaks))
 cols <- heat.colors(length(brks)-1)
 oldpar <- par(mfrow=c(nrowspar,tlength/nrowspar))
  for (i in 1:tlength) {
  image(x.seq0,y.seq0,f[[1]][,,i],asp=1,xlab="",ylab="",main=i,breaks=brks,col=cols)
  contour(x=x.seq0,y=y.seq0,z=pvalue[,,i],levels = c(1-alpha,1.5),col ='red',lwd='2',add=T)  #p-value contour
  #contour(x=x.seq0,y=y.seq0,z=f[[1]][,,i],levels = c(1-alpha,1.5),col ='blue',lwd='2',add=T)  #kde contour
#title(main="Excess risk surface (bandwidth=0.04)")
#plot(guichi2,axes=TRUE,add=TRUE)
#lines(rivers2)
#points(point)
  }
 par(oldpar)
   }    
if (tlength%%2==1) {
  tlength2<-tlength+1
  brks <- quantile(f[[1]], seq(0,1,breaks))
 cols <- heat.colors(length(brks)-1)
 oldpar <- par(mfrow=c(nrowspar,tlength2/nrowspar))   
  for (i in 1:tlength) {
 image(x.seq0,y.seq0,f[[1]][,,i],asp=1,xlab="",ylab="",main=i,breaks=brks,col=cols)
 contour(x=x.seq0,y=y.seq0,z=pvalue[,,i],levels = c(1-alpha,1.5),col ='red',lwd='2',add=T) #p-value contour
 #contour(x=x.seq0,y=y.seq0,z=f[[1]][,,i],levels = c(1-alpha,1.5),col ='blue',lwd='2',add=T) #kde contour
  }
  par(oldpar)
   }
 return(list(dens=f[[1]],pvalue=pvalue))
}