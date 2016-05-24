#########################################################################################################
# 5. Multivariate glyph plot of GWPCA loadings ##########################################################
#########################################################################################################


# This could do with some manipulation so that the size of the glyphs are controllable...
glyph.plot <- function(ld,loc,
	r1=50,
	add=FALSE,alpha=1,sep.contrasts=FALSE) {
	
	r <- max(max(loc[,1])-min(loc[,1]),max(loc[,2])-min(loc[,2]))/r1

  glyph.plot1 <- function(ld,loc,r=max(max(loc[,1])-min(loc[,1]),max(loc[,2])-min(loc[,2]))/r1,add=FALSE,alpha=1) {
		rowmax <- function(z) z[cbind(1:nrow(z),max.col(abs(z)))]
		ld <- sweep(ld,1,sign(rowmax(ld)),'*')
		ld.max <- max(ld)
		ld.min <- min(ld)
		n.col <- ncol(ld)
		n.row <- nrow(ld)
		angles <- (0:(n.col-1))*pi/n.col
		J <- 0+1i
		disp <- exp((pi/2-angles)*J)*r
		loc2 <- loc[,1] + loc[,2]*J
		ld.scaled <- 2*ld/(max(ld.max,-ld.min))
		if (!add) plot(loc,asp=1,type='n')
		points(loc,pch=16,cex=0.1,col='black')
		for (i in 1:n.row) {
			for (j in 1:n.col) {
				l.from <-  loc2[i]
				l.to   <-  loc2[i]+disp[j]*ld.scaled[i,j]
				col <- if (ld[i,j] > 0) {rgb(0,0,1,alpha)}
					   else {rgb(1,0,0,alpha)}
				lines(Re(c(l.from,l.to)),Im(c(l.from,l.to)),col=col) }}}

	glyph.plot2 <- function(ld,loc,r=max(max(loc[,1])-min(loc[,1]),max(loc[,2])-min(loc[,2]))/r1,add=FALSE,alpha=1) {
		rowmax <- function(z) z[cbind(1:nrow(z),max.col(abs(z)))]
		ld <- sweep(ld,1,sign(rowmax(ld)),'*')
		ld.max <- max(ld)
		ld.min <- min(ld)
		n.col <- ncol(ld)
		n.row <- nrow(ld)
		angles <- (0:(n.col-1))*2*pi/n.col
		J <- 0+1i
		disp <- exp((pi/2-angles)*J)*r
		loc2 <- loc[,1] + loc[,2]*J
		ld.scaled <- abs(ld)/(max(ld.max))
		if (!add) plot(loc,asp=1,type='n')
		points(loc,pch=16,cex=0.1,col='black')
		for (i in 1:n.row) {
			for (j in 1:n.col) {
				l.from <-  loc2[i]
				l.to   <-  loc2[i]+disp[j]*ld.scaled[i,j]
				col <- if (ld[i,j] > 0) {rgb(0,0,1,alpha)}
					   else {rgb(1,0,0,alpha)}
				lines(Re(c(l.from,l.to)),Im(c(l.from,l.to)),col=col) }}}

		if (sep.contrasts)
			{glyph.plot1(ld,loc,r,add,alpha)}
		else
			{glyph.plot2(ld,loc,r,add,alpha)}}


# To interact with the glyph map...
check.components <- function(ld,loc) {
	for (i in 1:ncol(ld)) cat(sprintf("%10s ",colnames(ld)[i]))
	cat('\n')
	repeat {
		this <- identify(loc,n=1,plot=FALSE)
		if (length(this)==0) break
		for (i in 1:ncol(ld)) cat(sprintf("%10.5f ",ld[this,i]))
		cat('\n')}
    }
    
    
########################################GW PCP plot
# GW parallel coordinate plot (PCP) to inspect an individual observation against its neighbouring observations
# The transparency of the neighbouring observation plot lines increases with distance...
# Can use this plot to inspect for a possible outlier (as identifed above)...
# Note there are many parameters (esp for 'wts' and 'tsc') that can be tweaked with this function
# & thus it needs to be made more general...
gw.pcplot <- function(data,vars,focus,bw,adaptive = FALSE, ylim=NULL,ylab="",fixtrans=FALSE, p=2, theta=0, longlat=F,dMat,...) 
{
	if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    loc<-coordinates(data)
  }
  else
     stop("Given data must be a Spatial*DataFrame")
  data <- as(data, "data.frame")
  dp.n<-nrow(data)
  i<-focus
  col.nm<-colnames(data)
  var.idx<-match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if(length(var.idx)==0) stop("Variables input doesn't match with data")
  x<-data[,var.idx]
  x<-as.matrix(x)
  m <- ncol(x)
  if (missing(dMat))
  {
    DM.given<-F
    if(dp.n <= 5000)
    {
      dMat <- gw.dist(dp.locat=loc, p=p, theta=theta, longlat=longlat)
      DM.given<-T
    }
  }
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
       stop("Dimensions of dMat are not correct")
  }
    
	if(DM.given)
	   dists <- dMat[,i]
  else
     dists <-gw.dist(dp.locat=loc, focus = i, p=p, theta=theta, longlat=longlat)
  if(adaptive)
  {
    rnk<- rank(dists, ties.method='first')
		bw<- dists[rnk==bw]  
  }
  #dists <- dists**2
	nbrlist <- which(dists < bw)
  dists <- dists**2
  wts <- (1 - dists/(bw*bw))^12
#  if(adaptive)
#    bw<-dists[bw]
  #nbrlist <- which(dists < bw*bw)
	xss <- scale(x)
	span <- 1:m
	tsc <- 25/length(nbrlist)
	if (is.null(ylim)) 
     ylim <- c(min(xss[nbrlist,]),max(xss[nbrlist,]))
	plot(span,xss[i,],type='l',ylim=ylim,
		xlim=c(0.5,m+0.5),col='red',lwd=6,axes=FALSE,xlab="",ylab=ylab,...)
	axis(1,at=1:m,labels=colnames(x),las=2,cex.axis=1.2)
	axis(2,at=seq(floor(ylim[1]),ceiling(ylim[2]),by=1),cex.axis=1.2)
	abline(v=1:m,col=grey(0.6))
	lines(c(1,m),c(0,0),col=grey(0.6))
	if (fixtrans) {
	    for (nbr in nbrlist) 
          lines(span,xss[nbr,],col=rgb(0,0,0,0.3), lwd=3) }
  else 
  {
      for (nbr in nbrlist)
          lines(span,xss[nbr,],col=rgb(0,0,0,tsc*wts[nbr]), lwd=3)  
  } 

}