stkde<-function(xlong,ylat,ztime,xgrids,ygrids,breaks=0.05,alpha=0.05,nrowspar=1,...) {
  #=======================================================
  #
  #  TITLE:     Spatio-Temporal Kernel Density Estimation with density contours
  #  FUNCTION:  stkde()
  #  AUTHOR:    Zhijie Zhang
  #  DATE:      01 JANUARY 2010
  #  CALLS:     np,graphics
  #  NEEDS:
  #  NOTES:
  # xlong-Projected planar coordinates of longitude
  # ylat- Projected planar coordinates of latitude
  # ztime-the integer variable,such as YEAR
  # xgrids,ygrids-number of grids to evaluate the density in the x and y directions
  # samkde<-stkde(xlong=d$x,ylat=d$y,ztime=d$tf,xgrids=10,ygrids=10,breaks=0.05,bwmethod="cv.ml")
  #  CHANGE HISTORY:
  #=======================================================
 # SET DEPENDENCIES
 require(np)
 #basic data manipulation
 year.seq <- sort(unique(ztime))
 tlength<-length(year.seq)
 #x.seq <- seq(0,1,length=xgrids)
 #y.seq <- seq(0,1,length=ygrids)
 x.seq <- seq(floor(min(xlong)),ceiling(max(xlong)),length=xgrids)
 y.seq <- seq(floor(min(ylat)),ceiling(max(ylat)),length=ygrids)
 #bandwidth selection and density estimation
 bw <- npudensbw(formula=~xlong+ylat+ordered(ztime),...)
 data.eval <- expand.grid(xlong=x.seq,ylat=y.seq,ztime=year.seq)
 fgrid<-npudens(bws=bw, newdata=data.eval)
 #bb<-data.frame(fgrid$eval,fgrid$dens)
 #fhat <- fitted(npudens(bws=bw, newdata=data.eval))
 fhat <-fgrid$dens
 f <- array(fhat,c(xgrids,ygrids,tlength))
 #mapping the results
 #brks <- quantile(f, seq(0,1,breaks))
 #cols <- heat.colors(length(brks)-1)
 #oldpar <- par(mfrow=c(1,tlength))
 #for (i in 1:tlength) image(x.seq,y.seq,f[,,i],asp=1,xlab="",ylab="",main=i,breaks=brks,col=cols)
 #par(oldpar)
 #Now we can map the KDE map and contours
if (tlength%%2==0)  {
 #mapping the results,even number
 brks <- quantile(f, seq(0,1,breaks))
 cols <- heat.colors(length(brks)-1)
 oldpar <- par(mfrow=c(nrowspar,tlength/nrowspar))  #change the number of rows and colunms for a map
  for (i in 1:tlength) {
  image(x.seq,y.seq,f[,,i],asp=1,xlab="",ylab="",main=i,breaks=brks,col=cols)
  contour(x=x.seq,y=y.seq,z=f[,,i],levels = c(1-alpha,1.5),col ='blue',lwd='2',add=T)  #kde contour
#probability=density*area,here there is problem for the results,maybe have been done this by that function in np package
#title(main="Excess risk surface (bandwidth=0.04)")
#plot(guichi2,axes=TRUE,add=TRUE)
#lines(rivers2)
#points(point)
  }
 par(oldpar)
   }    
if (tlength%%2==1) {
  tlength2<-tlength+1
  brks <- quantile(f, seq(0,1,breaks))
 cols <- heat.colors(length(brks)-1)
 oldpar <- par(mfrow=c(nrowspar,tlength2/nrowspar))  #change the number of rows and colunms for a map 
  for (i in 1:tlength) {
 image(x.seq,y.seq,f[,,i],asp=1,xlab="",ylab="",main=i,breaks=brks,col=cols)
 contour(x=x.seq,y=y.seq,z=f[,,i],levels = c(1-alpha,1.5),col ='blue',lwd='2',add=T) #kde contour
  }
  par(oldpar)
   }
 return(list(bw=bw$bw,dens=f))
}