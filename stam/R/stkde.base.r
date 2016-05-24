stkde.base<-function(xlong,ylat,ztime,xgrids,ygrids,breaks=0.05,...) {
  #=======================================================
  #
  #  TITLE:     Spatio-Temporal Kernel Density Estimation
  #  FUNCTION:  stkde.base()
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
 #return(list(bw=bw$bw,dens=f))
 #return(list(dens=f))
}