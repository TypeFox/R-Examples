plotsML.gb2 <- function(z, shape1, scale, shape2, shape3, w=rep(1,length(z))){
  
  d     <- data.frame(inc=z,w=w)
	d     <- d[!is.na(d$inc),]
  index <- order(d$inc)
  d     <- d[index,]
  inc   <- d$inc
  w     <- d$w
  cdf   <- (cumsum(w) - w/2)/sum(w)

# Truncate at 0

  d     <- data.frame(inc=inc, w=w, cdf=cdf ) 
  x     <- d$inc[d$inc > 0]
  cdf   <- d$cdf[d$inc > 0]
  w     <- d$w[d$inc > 0]
  
  par(mfrow=c(2,1))

  limi=0
  kernel="epanechnikov"
  shift <- min(cdf[x>=limi])

# Computation of the indicators in the original scale

  fshift <- function(p, shift) (p-shift)/(1-shift)

# Inverse transformation
  ifshift <- function(px, shift) {shift + px*(1-shift)}

  #limsup <- qgb2(0.999,shape1,scale,shape2,shape3)
  limsup <- max(x[cdf<= 0.999])+limi
  #liminf  <- qgb2(0.001,shape1,scale,shape2,shape3)
  liminf <- 0
  xlim <- c(liminf,limsup)+limi	

# Cumulative distribution plot
	plot(cdf,ifshift(pgb2(x,shape1,scale,shape2,shape3),shift), type="l", col="red", main="Cumulative Distribution plot",
	xlab="empirical distribution", ylab="GB2 distribution")
	abline(0,1,lty=2)

# Density plot
	weights <- w/sum(w)	
	dens <- density(x, weights=weights, kernel=kernel)

# GB2 mode
	modeGB2 <- scale*((shape1*shape2-1)/(shape1*shape3+1))^(1/shape1)
	ysup <- max(dgb2(modeGB2, shape1, scale, shape2, shape3), max(dens$y))
	leg.txt <- c(paste("Kernel", "   ", sep=""), "GB2 density")
  xx1 <- qgb2(0.9, shape1, scale, shape2, shape3)
	xx2 <- qgb2(0.6, shape1, scale, shape2, shape3)
	yy <- dgb2(xx2, shape1, scale, shape2, shape3)
	plot(dens,main=paste("Density plot")#, sub=paste("kernel=",kernel)
	,xlim=xlim,ylim=c(0,ysup))
	curve(dgb2(x,shape1,scale,shape2,shape3),from=liminf,to= limsup,n=2000, col="red", add = TRUE)
	legend(list(x=xx1,y=yy), legend = leg.txt, col=c("black", "red"), lty=1, bty="n")
}

saveplot <- function(name, pathout){
 	local({
 		dev.set (2)
 		nameplot <- paste(pathout,name,".pdf",sep="")
 		dev.print (device=pdf, file=nameplot);
 	      })
}
