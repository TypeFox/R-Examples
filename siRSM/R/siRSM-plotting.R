#################################################################################
# surf.stats: a support function, computing stationary points, principal axes
# for response surface a la Edwards (2002) 
#################################################################################
surface.stats.main <- function(b,xlim,ylim)
{
# Written by Huan Cheng, 2014
# Modified by Mu Zhu, 2014
# b = coefficients of quadratic response surface
# xlim, ylim = the boundary on the floor in the 3-d graph 
  b0 <- b[1] ; b1 <- b[2]; b2 <- b[3]
  b3 <- b[4] ; b4 <- b[5]; b5 <- b[6]
  
  ### calculate the stationary point ... a la Edwards (2002) 
  if (4*b3*b5-b4^2 == 0 ) stop("denominator is zero")
  stnx <- as.vector((b2*b4-2*b1*b5)/(4*b3*b5-b4^2))
  stny <- as.vector((b1*b4-2*b2*b3)/(4*b3*b5-b4^2))
  
  ### calculate the principal axes ... a la Edwards (2002)
  p11 <- as.vector((b5-b3+sqrt((b5-b3)^2+b4^2))/b4)
  p10 <- as.vector(stny-stnx*p11)
  p21 <- as.vector((b5-b3-sqrt((b5-b3)^2+b4^2))/b4)
  p20 <- as.vector(stny-stnx*p21)
 
  ### calculate ax, ax2 along congruence and incongruence line ... a la Edwards (2002)
  ax.congr <- as.vector(b1+b2); ax2.congr <- as.vector(b3+b4+b5)
  ax.incongr <- as.vector(b1-b2); ax2.incongr <- as.vector(b3-b4+b5)
  
 if (!missing(xlim) && !missing(ylim)){
  ### compute where the axes intersect boundaries of plot
  pxl <- xlim[1]
  pyl <- p10+p11*pxl
  if (pyl < ylim[1]) {
    pyl <- ylim[1]; pxl <- (pyl-p10)/p11
  }
  if (pyl > ylim[2]) {
    pyl <- ylim[2]; pxl <- (pyl-p10)/p11
  }
  pxh <- xlim[2]
  pyh <- p10+p11*pxh
  if (pyh > ylim[2]){
    pyh <- ylim[2]; pxh <- (pyh-p10)/p11
  }  
  if (pyh < ylim[1]){
    pyh <- ylim[1]; pxh <- (pyh-p10)/p11
  }  
  pl <- as.vector(c(pxl,pyl))
  ph <- as.vector(c(pxh,pyh))
   
  sxl<- xlim[1]
  syl <- p20+p21*sxl
  if(syl > ylim[2]){
    syl <- ylim[2]; sxl <- (syl-p20)/p21
  }
  if(syl < ylim[1]){
    syl <- ylim[1]; sxl <- (syl-p20)/p21
  }
  sxh <- xlim[2]
  syh <- p20+p21*sxh
  if (syh < ylim[1]){
    syh <- ylim[1]; sxh <- (syh-p20)/p21
  }  
  if (syh > ylim[2]){
    syh <- ylim[2]; sxh <- (syh-p20)/p21
  }  
  sl <- as.vector(c(sxl,syl))
  sh <- as.vector(c(sxh,syh))
 }
 else{
  pl<-ph<-sl<-sh<-NULL
 } 
  return(list(u0=stnx, v0=stny,
              p10=p10, p11=p11,
              p20=p20, p21=p21,
			  ax.congr=ax.congr, ax2.congr=ax2.congr,
			  ax.incongr=ax.incongr, ax2.incongr=ax2.incongr,
			  pl=pl,ph=ph,
              sl=sl,sh=sh))
}

surface.stats <- function(obj)
{
  if (obj$type=='interaction.only'){
    b=c(obj$coef[1],obj$coef[2],obj$coef[3],0,obj$coef[4],0)
  }
  else{
    b=obj$coef
  }
  tmp=surface.stats.main(b)
  return(data.frame(u0=round(tmp$u0,4), v0=round(tmp$v0,4),
              p10=round(tmp$p10,4), p11=round(tmp$p11,4),
              p20=round(tmp$p20,4), p21=round(tmp$p21,4),
			  ax.congr=round(tmp$ax.congr,4), ax2.congr=round(tmp$ax2.congr,4),
			  ax.incongr=round(tmp$ax.incongr,4), ax2.incongr=round(tmp$ax2.incongr,4)))
}

#################################################################################
# plot.siRSM: creates 3D perspective plot of the siRSM model, uses 'rsm' package
#################################################################################
draw.full.quadratic <- function(x, xname=NULL, yname=NULL, zname=NULL, center='zero', debug=FALSE)
{
# Written by Huan Cheng, 2014
# Modified by Mu Zhu, 2014
  dxyz <- data.frame(u=x$u, v=x$v, y=x$y)
  lm.full <- lm(y~u+v+I(u^2)+I(u*v)+I(v^2),data=dxyz)
  max.x <- max(x$u); min.x <- min(x$u) 
  max.y <- max(x$v); min.y <- min(x$v)

  ### extend plotting space so that stationary point is included
  tmp <- surface.stats.main(coef(lm.full))
  st=c(tmp$u0,tmp$v0)
  if (st[1] < min.x) { min.x = st[1] - (max.x-st[1]) }
  if (st[1] > max.x) { max.x = st[1] + (st[1]-min.x) }
  if (st[2] < min.y) { min.y = st[2] - (max.y-st[2]) }
  if (st[2] > max.y) { max.y = st[2] + (st[2]-max.y) }  
  r.x=max.x-min.x; r.y=max.y-min.y
  min.x = min.x - r.x/10; max.x = max.x + r.x/10
  min.y = min.y - r.y/10; max.y = max.y + r.y/10
  min.x<-min.y<-min(min.x,min.y)
  max.x<-max.y<-max(max.x,max.y)
  if (center=='zero') {
 	R=max(abs(min.x),abs(max.x))
	min.x<-min.y<-(-R)
	max.x<-max.y<-R
  }  
  
  ### Draw the surface plot 
  res <- persp(lm.full,
               v~u, zlab='', col='yellow',
			   bounds=list(u=c(min.x,max.x),v=c(min.y,max.y)),
			   contours=list(z='bottom'), cex.lab=1.5) 
  if (debug) print(res)

  ### find stationary point and two principle axes ... a la Edwards (2002)
  st <- surface.stats.main(coef(lm.full),c(min.x,max.x),c(min.y,max.y))
  if (debug) print(st)
  xx1 <- st$pl
  xx2 <- st$sl
  yy1 <- st$ph
  yy2 <- st$sh
  p1 <- st$u0
  p2 <- st$v0
  minv <- min(res[[1]]$zlim)
  res <- res[[1]]$transf

  # data points going into the fitting of response surface
  points(trans3d(x$u,x$v,z=minv,pmat=res),cex=0.5,col='grey')
  
  # 45 degree lines 
  xlxh <- c(min.x,max.x)
  ylyh <- c(min.y,max.y)
  yhyl <- c(max.y,min.y) 
  lines(trans3d(xlxh,ylyh,z=minv,pmat=res),lty=3)
  lines(trans3d(xlxh,yhyl,z=minv,pmat=res),lty=3)
  
  # x- and y-axes
  lines(trans3d(c(0,0),c(min.y,max.y),z=minv,pmat=res),lty=3)
  lines(trans3d(c(min.x,max.x),c(0,0),z=minv,pmat=res),lty=3)
  
  # principal axes
  lines(trans3d(c(xx1[1],yy1[1]),c(xx1[2],yy1[2]),z=minv,pmat=res),lty=1,col="red",lwd=2)
  lines(trans3d(c(xx2[1],yy2[1]),c(xx2[2],yy2[2]),z=minv,pmat=res),lty=1,col="blue",lwd=2)
  
  # stationary point
  points(trans3d(p1,p2,z=minv,pmat=res),pch=19,cex=2)
  
  # title/label
  if (!is.null(xname) && !is.null(yname) && !is.null(zname)){ 
    title(main=paste('y=f(u,v)\n','u = ',xname,', v = ',yname,', y = ',zname,sep=''),cex.main=1.5)
  }
  else{
    title(main='y=f(u,v)',cex.main=1.5)
  }
  
  subtitle='Solid Red = 1st Principal Axis; Solid Blue = 2nd Principal Axis\nDashed: u=0, v=0, u+v=0, and u-v=0'
  title(sub=subtitle,cex=1.25)
}

draw.interaction.only <- function(x, xname=NULL, yname=NULL, zname=NULL)
{
  dxyz <- data.frame(u=x$u, v=x$v, y=x$y)
  lm.int <- lm(y~u+v+I(u*v),data=dxyz)
  u.lo=mean(x$u)-sd(x$u)
  u.hi=mean(x$u)+sd(x$u)
  v.lo=mean(x$v)-sd(x$v)
  v.hi=mean(x$v)+sd(x$v)
  # four prototypical observations
  # (lo, lo)
  # (lo, hi)
  # (hi, lo)
  # (hi, hi)
  newdata=cbind(
   c(u.lo,u.lo,u.hi,u.hi),
   c(v.lo,v.hi,v.lo,v.hi))
  newdata=as.data.frame(newdata)
  dimnames(newdata)[[2]]<-c('u','v')
  newy=predict(lm.int,newdata=newdata)

  rx=max(newdata[,1])-min(newdata[,1])
  x.range=c(min(newdata[,1])-rx/10,max(newdata[,1])+2*rx/10)
  ry=max(newy)-min(newy)
  y.range=c(min(newy)-ry/10,max(newy)+ry/10)  
  plot(newdata[,1],newy,axes=F,
    xlim=x.range,ylim=y.range,
    type='n',xlab='',ylab='y=f(u,v)',cex.lab=1.35)
  points(c(u.lo,u.hi),c(newy[1],newy[3]),pch=15,cex=2,col='blue')
  points(c(u.lo,u.hi),c(newy[2],newy[4]),pch=19,cex=2,col='red')
  segments(u.lo,newy[1],u.hi,newy[3],lty=2,lwd=2,col='blue')
  segments(u.lo,newy[2],u.hi,newy[4],lty=1,lwd=2,col='red')
  text(u.hi+rx/10,newy[3],'v=LO',cex=1.35)
  text(u.hi+rx/10,newy[4],'v=HI',cex=1.35)
  axis(1,at=c(u.lo,u.hi),labels=c('u=LO','u=HI'),cex.axis=1.35)
  axis(2,cex.axis=1.35)
  box()

  if (!is.null(xname) && !is.null(yname) && !is.null(zname)){ 
    title(main=paste('u = ',xname,', v = ',yname,', y = ',zname,sep=''),cex.main=1.5)
  }
}

plot.siRSM.default <- function(x, ...)
{
  if (x$type=='interaction.only') {
    draw.interaction.only(x,...)
  }
  else{
    draw.full.quadratic(x,...)
  }  
}

plot.siRSM <- function(x, ...) 
UseMethod("plot.siRSM")

