
# Fibroblast example as in Section 4.2 of Advances in dendrogram seriation for application to visualization, D. Earle and C. Hurley

library(DendSer)




if (!("Iyer517" %in% rownames(installed.packages()))) install.packages("Iyer517")

library(Iyer517)
data(Iyer517)
fib <- log(exprs(Iyer517)[,1:12],2)
#This is what Eisen et al. (1998) do and this also gives
#the data that Tien et al. (2008) use.



colnames(fib) <- c("0","15m","30m","1hr","2hr","3hr","4hr","8hr","12hr","16hr","20hr","24hr")



dm <- 1- cor(t(fib))
d <- as.dist(dm)
h <- hclust(d,method="ward")

ob <- dser(h,dm,cost=costBAR)
ob<-rev(ob)


od <- dser(h, ser_weight=rowMeans(fib),cost=costLS)  
od<-rev(od)

#---------------------------
# plots of corr matrix



mcol <-sequential_hcl(50)


dev.new(height=3, width=3);par(mar=c(1,1,1,1))
plotAsColor(dm,ann=F,order=od,useRaster=T,col=mcol) # Fig 11a


dev.new(height=3, width=3);par(mar=c(1,1,1,1))
plotAsColor(dm,ann=F,order=ob,useRaster=T,col=mcol) # Fig 11b

#---------------------------


# plots of data 

breaks <- c(min(fib)-1,seq(-3,3,length.out=30),max(fib)+1)
hcol <-diverge_hcl(length(breaks)-1)


# dev.new(width=2.35)
# par(mar=c(3,1,1,1))
# plotAsColor(fib[od,],col=hcol,breaks=breaks,useRaster=T)
# axis(1, at=c(1:12), las=2, labels=colnames(fib))


dev.new(width=2.35) # Fig 12a
par(mar=c(3,1,1,1))
plotAsColor(fib[ob,],col=hcol,breaks=breaks,useRaster=T)
axis(1, at=c(1:12), las=2, labels=colnames(fib)) 



nc<-10
ccols <- rainbow(nc)


clush <- cutree(h, nc)
names(clush)<-NULL
zz<- tapply(1:nrow(fib),clush[ob],min )
clus<- rank(zz)[clush]   # renumber the clusters, consistent with order in ob



# cluster legend for data ordered by ob

dev.new(width=.5)  # Fig 12a- legend
par(mar=c(3,1,1,.5))
plotAsColor(as.matrix(clus),order.row=ob,col=ccols[1:nc],useRaster=T)
clusname<-LETTERS[1:nc]
z <- nrow(fib) -tapply(order(ob),clus,mean)
axis(2,at=z,labels=clusname,las=2,tick=F,line=-.6,cex.axis=.7)


# cluster legend for data ordered by od

dev.new(width=.5)
par(mar=c(3,1,1,.5))

plotAsColor(as.matrix(clus),order.row=od,col=ccols[1:nc],useRaster=T)
clusname<-LETTERS[1:nc]
z <- nrow(fib) -tapply(order(od),clus,mean)
axis(2,at=z,labels=clusname,las=2,tick=F,line=-.6,cex.axis=.7)



#---------------------------
fibm<-aggregate(fib, list(clus),mean)

# Fig 12b
dev.new(width=3.5)
  m<- rbind(matrix(1:10, 5, 2, byrow = F),c(11,12))
  layout(m,widths = c(1,1),heights=c(1,1,1,1,1,.2))
   par(mar=c(1,3,0,.3))  
   par(mgp=c(3, 1, 0))
    for (i in 1:nc) {
  	plot(1:12,fibm[i,-1], col=ccols[i],ylim=c(-1.5,2.1),axes=F,pch=20,ann=F)
   	lines(1:12,fibm[i,-1], col=ccols[i])
  	axis(2, at=c(-1.5,1.5),las=2,cex.axis=.7)
  	 mtext(clusname[i],2,line=.75,las=2,cex=.7)
  	abline(h=0,col="grey50")
  }
   plot(1:12,fibm[i,-1], col=ccols[i],ylim=c(-1.5,2.1),axes=F,main="",type="n",ann=F)
  par("mgp" = c(3,-.5,1),"tcl"= .2)
 axis(3,at=1:12,labels=colnames(fib),las=2,cex.axis=.7)
 plot(1:12,fibm[i,-1], col=ccols[i],ylim=c(-1.5,2.1),axes=F,main="",type="n",ann=F)
  par("mgp" = c(3,-.5,1),"tcl"= .2)
 axis(3,at=1:12,labels=colnames(fib),las=2,cex.axis=.7)




# clusters on mds plot
mds2 <- cmdscale(d, k=2)
dev.new()
plot(mds2, pch=clusname[clus], col=ccols[clus], asp=1, main="2-D Classical MDS solution", 
      xlab="Dimension 1", ylab="Dimension 2")



