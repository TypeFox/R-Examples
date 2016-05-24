#get the cross points of lines a=slop, b=intercept
linescrosspoint<-function(a, b){
  b1=b[-NROW(b)]
  b2=b[-1]
  a1=a[-NROW(a)]
  a2=a[-1]
  
  x1=-(b1-b2)/(a1-a2)
  y1=-(a2*b1-a1*b2)/(a1-a2)
  na.exclude(data.frame(x=x1, y=y1))
}

#get range quantile for scatter of x, y
xy.quantile<-function(x, y, alpha=0.05, half.nknots=90) { 
  probs=range(c(alpha/2, 1-alpha/2))
  anglelist=seq(-pi/2, pi/2-pi*(0.5/half.nknots), length.out=half.nknots+1)
  res=c(1)[0]
  off.points=rep(0, NROW(x))
  for (angle in anglelist) {
    if (angle%%(pi/2)!=0) {
      a=tan(angle)
      b=-1
      d=-1*(a*x+b*y)/sqrt(a^2+b^2)
      crt=quantile(x=d, probs=probs)
      idx=which((d<crt[1])|(d>crt[2]), arr.ind=T)
      off.points[idx]=1
      res=rbind(res, c(a, crt*sqrt(1+a^2), angle, 1))
    }
  }
  
  
  lines=data.frame(b=c(res[,1],res[,1]), incpt=c(res[,2],res[,3]), angle=c(res[,4], res[,4]+pi)%%(2*pi), c=c(res[,5], res[,5]))
  
  p=linescrosspoint(a=c(lines[,1], lines[1:2,1]), b=c(lines[,2], lines[1:2,2]))
  list(polygon=p, data=data.frame(x=x, y=y, pip=1-off.points), alpha.realized=sum(off.points)/NROW(x))
}

plot.distfree.cr<-function(x, show.points=T, ...) {
obj=x;
  xlim=range(c(obj$data[,1], obj$polygon$x))
  ylim=range(c(obj$data[,2], obj$polygon$y))
  plot(NA, type="n", xlim=xlim, ylim=ylim, xlab=obj$xlab, ylab=obj$ylab)
  if (show.points) {
    points(obj$data, col=obj$col[2:3][2-obj$data$pip])
  }
  polygon(obj$polygon, border=obj$col[1], lwd=2.0)
  polygon(obj$polygon.smooth1, border=obj$col[1], lwd=2.0, lty="dotted")
  polygon(obj$polygon.smooth2, border=obj$col[1], lwd=2.0,  lty="dashed")
  legend("topright", legend=c("Apex connection", "Smooth1", "Smooth2"), col=obj$col[1], lty=c("solid", "dotted", "dashed"), bty="n")
}

distfree.cr<-function(x, y, alpha=0.05, alpha.min.diff=alpha/10, nknots=40,
                      xlab = deparse(substitute(x)), 
                      ylab = deparse(substitute(y)),
                      col=c("red", "black", "gray"), draw=T) {
  if (missing(y)) {
    if (NCOL(x) == 2) {
      if (missing(xlab)) xlab <- colnames(x)[1]
      if (missing(ylab)) ylab <- colnames(x)[2]
      y <- x[, 2]
      x <- x[, 1]
    } else stop("x and y must be vectors, or x must be a 2 column matrix")
  } else if (!(is.vector(x) && is.vector(y) && length(x) == length(y))) stop("x and y must be vectors of the same length")
  col=rep(col, 3)
  
  est.alpha=alpha/nknots
  diff=0.5
  while (diff>alpha.min.diff) {
    xy.quant=xy.quantile(x=x, y=y, half.nknots=round(nknots/2), alpha=est.alpha)
    diff=abs(alpha-xy.quant$alpha.realized)
    est.alpha=est.alpha*alpha/xy.quant$alpha.realized
  }
  xy.quant$alpha=alpha
  xy.quant$xlab=xlab
  xy.quant$ylab=ylab
  xy.quant$col=col[1:3]
  
  #   #smooth the regions
  #   for (k in 2:(NROW(xy.quant$polygon)-1)) {
  #     if (sum((xy.quant$polygon[k-1,]-xy.quant$polygon[k,])^2)>sum((xy.quant$polygon[k-1,]-xy.quant$polygon[k+1,])^2)) {
  #       xy.quant$polygon[c(k, k+1),]=xy.quant$polygon[c(k+1, k),]  
  #     }
  #   }  
  
  xy.quant$polygon.smooth1=xy.quant$polygon;
  for (k in 2:(NROW(xy.quant$polygon.smooth1)-1)) {
    if (sum((xy.quant$polygon.smooth1[k-1,]-xy.quant$polygon.smooth1[k,])^2)>sum((xy.quant$polygon.smooth1[k-1,]-xy.quant$polygon.smooth1[k+1,])^2)) {
      xy.quant$polygon.smooth1[c(k, k+1),]=xy.quant$polygon.smooth1[c(k+1, k),]  
    }
  }
  
  

  idx=which(xy.quant$data$pip==1)
  tmp=xy.quant$data[,1:2][idx,];
  idx=border.points.idx(tmp)
  
  xy.quant$polygon.smooth2=tmp[idx,]
  
  class(xy.quant)<-"distfree.cr"
  if (draw) plot(xy.quant)
  xy.quant
}


#obtain angle between (0,0) and (x, y)
xytoangle<-function(x, y) {
  d=sqrt(x^2+y^2)
  ang.cos=acos(x/d)%%(2*pi)
  idx=which(y<0,arr.ind=T)
  ang.cos[idx]=2*pi-ang.cos[idx]
  ang.cos
}

#get index of outer points in x=data.frame(row, col)
border.points.idx<-function(x) {
  x=data.frame(row=x[,1], col=x[,2])
  
  k=order(x$row)[1]
  bord.points=k;
  ang0=pi/2;
  
  while(NROW(bord.points)==NROW(unique(bord.points))) { 
    ang=xytoangle(x=x$row-x$row[k], y=x$col-x$col[k])
    ang_diff=ang-ang0
    idx=which(ang_diff>pi,arr.ind=T) 
    ang_diff=ang_diff%%(2*pi)
    
    minang=min(ang_diff,  na.rm = T)
    k=which(minang==(ang_diff), arr.ind=T)
    bord.points=c(bord.points, k)
    ang0=ang[k]%%(2*pi)
  }
  bord.points
}

# dat=pt$dat[which(pt$dat[,3]==1),][,1:2]
# idx=border.points.idx(dat)
# polygon(dat[idx,][,1:2])

