
hexpie = function(x, y = NULL, z = NULL, n = 24,  shape = "hex",
 p.rule = "radial", decr.by.rank = NULL,  freq.trans = I, alpha.freq = FALSE, col = "hcl",
 col.opt = list(), show.hex = TRUE, random = NULL, xlim = range(x), ylim = range(y),
 label.opt = list(), vp = NULL){

#c("angle","radius","alpha", "sample")
#require(hexbin)

x.name <- names(x)
y.name <- names(y)
if(is.null(x.name)) x.name<-"x"
if(is.null(y.name)) y.name<-"y"

if( "bg" %in% names(col.opt) ){
	bg <- col.opt$bg
}else{
	bg <- NA# hsv(0,0,0, alpha=0.02)#hcl(240,10,85)
}

if( "bgs" %in% names(col.opt) ){
	bgs <- col.opt$bgs
}else{
	bgs <- hsv(0,0,0, alpha=0.05)
}

if( "col.axis" %in% names(col.opt) ){
    col.axis <- col.opt$col.axis
}else{
    col.axis <- 1
if(!is.na(bg)){
    if(bg == 1 | bg == "black") col.axis <- "white"
}
}

if( "alpha.r" %in% names(col.opt) ){
	stopifnot( col.opt$alpha.r <= 1 & col.opt$alpha.r > 0 )
	alpha.r <- -log(col.opt$alpha.r)
}else{
	alpha.r <- -log(0.2)
}

if( "line.col" %in% names(col.opt) ){
	line.col <- col.opt$line.col
}else{
	line.col <- NA
}
if( "line.col.hex" %in% names(col.opt) ){
	line.col.hex <- col.opt$line.col.hex
}else{
	line.col.hex <- "darkgrey"
}
if( "alpha.hex" %in% names(col.opt) ){
	alpha.hex <- col.opt$alpha.hex
}else{
	alpha.hex <- FALSE
}

dbin <- hexbin(x,y, xbins = n, IDs = TRUE, xbnds = xlim, ybnds = ylim)

xy <- hcell2xy(dbin)

z <- as.integer(as.factor(z))

#xmarks <- unique(sort(xy$x))
#ymarks <- unique(sort(xy$y))


xrad <- diff(dbin@xbnds)/dbin@xbins/2
yrad <- diff(ylim)/n/sqrt(3)

xrng <- diff(xlim)
yrng <- diff(ylim)

xrng2 <- diff(range(xy$x))+2*xrad
yrng2 <- diff(range(xy$y))+2*yrad

xmarks <- seq(range(xy$x)[1]-xrad, range(xy$x)[2],2*xrad) #(xrng2-xrad)/(dbin@dimen[2]))
ymarks <- seq(range(xy$y)[1]-yrad, range(xy$y)[2]+yrad/2, 3*yrad/2)#yrng2/(dbin@dimen[1]-2))

xy <- cbind((xy$x-min(xy$x)+xrad)/xrng2,(xy$y-min(xy$y)+yrad)/yrng2)




data <- data.frame(cbind(dbin@cID,z))
data <- melt(data,1)
#data <- reshape2:::cast(data, V1 ~ value, fun.aggregate = length)
data <- acast(data, V1 ~ value, fun.aggregate = length)

if(!is.null(vp)){
	pushViewport(vp)	
}else{
	grid.newpage()	
}
grid.rect(gp=gpar(fill=bg))
vp.base <- viewport(0.54,0.46, width = 0.90, height = 0.90)
pushViewport(vp.base)

# hexagon
 angles = seq(30,330,60)/180*pi

#corners <-  cbind( cos(angles), sin(angles))/n/sqrt(3)
corners <-  cbind( cos(angles)*xrad/xrng2*2/sqrt(3), sin(angles)*yrad/yrng2)

if(show.hex){
	if(alpha.hex){
		#nix <- apply(cbind(xy,dbin@count/max(dbin@count)),1,function(s){
		#grid.polygon(x = (s[1]+corners[,1]), y = (s[2]+corners[,2]), gp = gpar(fill= alpha("black",sqrt(s[3]/2)), col = line.col.hex))})
	grid.polygon(rep(xy[,1],each=6)+corners[,1],
	rep(xy[,2],each=6)+corners[,2],
	id = rep(1:nrow(xy),each=6),
	gp=gpar(fill = alpha("black",sqrt(rep(xy[,31],each=6)/2)), col = line.col.hex))#'CHECK'

}else{
	grid.polygon(rep(xy[,1],each=6)+corners[,1],
	rep(xy[,2],each=6)+corners[,2],
	id = rep(1:nrow(xy),each=6),
	gp=gpar(fill = bgs, col = line.col.hex))
	
#nix <- apply(xy,1,function(s){
#	grid.polygon(x = (s[1]+corners[,1]), y = (s[2]+corners[,2]), gp = gpar(fill=bgs, col = line.col.hex))
	
#})	
}
}

angular <- any(c("angle","ang","a","angular","angles") %in% p.rule)
radial <- any(c("radius","rad","r","radial") %in% p.rule)


stopifnot(xor(angular, radial))

# prepare sizes
if(!is.null(freq.trans)){
	base.rad <- freq.trans(dbin@count)
	base.rad <- base.rad/max(base.rad)
	
}else{
	#base.rad<- rep(1,nrow(xy))
	base.rad <- dbin@count/max(dbin@count)
}
	 data <- data/dbin@count #data[,-1]/dbin@count #data$data
if(!is.null(random)){
	random <- as.integer(random)
	data <- t(apply(data,1,function(s){
		rmultinom(1,random,prob=s)	
		
	}))/random
	
}
DM <- cbind(xy,data, base.rad)

# instead: data == probs, base.rad == ntv+3 wts, xy == centers 



ntv <- ncol(DM)-3
colv <- getcolors(ntv,col,col.opt)

if(!is.null(decr.by.rank)){
	ofun <- function(s) order(s, decreasing = decr.by.rank)	
}else{
	ofun <- seq_along	
}


if(shape %in% c("c", "pie", "circular", "piechart") ){
	ngon.angles = seq(0,360,1/4)/180*pi
	ngon.corners <-  cbind( cos(ngon.angles)*xrad/xrng2, sin(ngon.angles)*xrad/xrng2)	
}



if( angular ){
		
		if(shape %in% c("hex", "hexagonal", "hexagon", "h") ){
			cat("Not implemented. Using radial shapes.")
			shape <- "pie"
		}
	
	if(shape %in% c("c", "pie", "circular", "piechart") ){
		
	
		apply(DM,1, function(s){
			xc <- s[1]
			yc <- s[2]
			probs <- s[-c(1,2,ntv+3)]
			rad <- ifelse(is.null(freq.trans),1,sqrt(s[ntv+3]))
		if(alpha.freq){
			colv2 <- alpha(colv,exp(alpha.r*(s[ntv+3]-1)))	
		}else{
			colv2 <- colv	
		}
			ord <- ofun(probs)
			cum.probs <-  rev(cumsum(probs[ord])) 
			mapply(function(s1,s2, s3){
			if(s3 > s1){
				cc = seq( s1, s3, (s3-s1)/ceiling( (s3-s1)*360 ))
				if(abs(s3-s1)==1){
					grid.polygon(x = xc + c(cos(cc*2*pi))*rad*xrad/xrng2, y = yc+c(sin(cc*2*pi))*rad*xrad/xrng2, gp = gpar(fill=s2, col = line.col))
				}else{
					grid.polygon(x = xc + c(0,cos(cc*2*pi))*rad*xrad/xrng2, y = yc+c(0,sin(cc*2*pi))*rad*xrad/xrng2, gp = gpar(fill=s2, col = line.col))
				}}}, s3 = cum.probs, s2 = colv2[rev(ord)], s1 = c(cum.probs[-1],0))
			return(invisible(TRUE))
		})	
	
	}
	
}




if( radial ){
	if(shape %in% c("hex", "hexagonal", "hexagon", "h") ){

		cum.probs <- #as.vector(
			apply(data,1,function(z){
				ord <- ofun(z)
				return( sqrt( rev(cumsum(z[ord])) ) )
			})#)
		
	
		cum.probs <- rep(cum.probs,each=6*2)
		
		if(alpha.freq){
			colv2 <- as.vector(
				sapply(base.rad,function(z){
					return( alpha(colv,exp(alpha.r*(z-1))) )
				}))
		}else{
			colv2 <- rep(colv,nrow(xy))
		}
		ord <- as.vector(
		apply(data,1,function(z){
				ord <- ofun(z)
				return( rev(ord) )
		}))
		ord <- ord + rep( (0:(nrow(xy)-1)),each=ntv)*ntv

		colv2 <- colv2[ord]
		
		colv2 <- rep(colv2, each = 2)
		colv2[seq(1,length(colv2),2)] <- "white"
		
		colv2 <- rep(colv2, 6)
	
		if(is.null(freq.trans)){
			rad <- 1
		}else{
			rad <- rep(sqrt(base.rad),each=ntv*6*2)
		} 
		
	
	
		grid.polygon(x = rep(xy[,1],each=2*6*ntv) + rep(c(corners[,1],corners[,1]),ntv)*rad*cum.probs,
			 y = rep(xy[,2],each=2*6*ntv)+rep(c(corners[,2],corners[,2]),ntv)*rad*cum.probs, 
			 gp = gpar(fill=colv2, col = line.col), id = rep(1:(nrow(xy)*ntv*2),each=6))
				

	}
	
	
	if(shape %in% c("c", "pie", "circular", "piechart") ){

		cum.probs <- #as.vector(
			apply(data,1,function(z){
				ord <- ofun(z)
				return( sqrt( rev(cumsum(z[ord])) ) )
			})#)
		
	
		cum.probs <- rep(cum.probs,each=2)
		
		if(alpha.freq){
			colv2 <- as.vector(
				sapply(base.rad,function(z){
					return( alpha(colv,exp(alpha.r*(z-1))) )
				}))
		}else{
			colv2 <- rep(colv,nrow(xy))
		}
		ord <- as.vector(
		apply(data,1,function(z){
				ord <- ofun(z)
				return( rev(ord) )
		}))
		ord <- ord + rep( (0:(nrow(xy)-1)),each=ntv)*ntv

		colv2 <- colv2[ord]
		
		colv2 <- rep(colv2, each = 2)
		colv2[seq(1,length(colv2),2)] <- "white"
		
		
	
		if(is.null(freq.trans)){
			rad <- 1
		}else{
			rad <- rep(sqrt(base.rad),each=ntv*2)
		} 
		
	
	
		grid.circle(x = rep(xy[,1],each=2*ntv) ,
			 y = rep(xy[,2],each=2*ntv), r = rad*cum.probs*xrad/xrng2, 
			 gp = gpar(fill=colv2, col = line.col))
				

	}
}

	 if( "cex" %in% names(label.opt) ){
		 label.opt$cex.axis <- label.opt$cex
		 label.opt$cex.axis <- label.opt$cex
	 }
	 if( "cex.axis" %in% names(label.opt) ){
		 cex.axis <- label.opt$cex.axis
	 }else{
		 cex.axis <- 1
	 }
	 if( "cex.lab" %in% names(label.opt) ){
		cex.lab <- label.opt$cex.lab
	 }else{
		 cex.lab <- 1
	 }
	 
popViewport()
vpx <- viewport(0.54,0.93,width=0.9,height = 0.05)
pushViewport(vpx)
my.grid.axis(x0=0,y0=0.1,len=1-xrad/xrng2,ticks=signif(xmarks,5), rot=0, keep = min(7,n-2), ltm=1/20 ,col.axis = col.axis, lab.cex=cex.axis, trot = 300 )
popViewport()
vpy <- viewport(0.07,0.46,width=0.05,height = 0.9)
pushViewport(vpy)
my.grid.axis(x0=0.9,y0=0,len=1-yrad/yrng2/2,ticks=signif(ymarks,5), rot=90, keep = min(7,n-2),ltm=1/20,col.axis = col.axis, lab.cex=cex.axis, trot = 340)			

popViewport()

vpyl <- viewport(0.02,0.5,width=0.05,height = 0.9)
pushViewport(vpyl)
grid.text(label=y.name,0.5,0.5,gp=gpar(fontsize = 14*cex.lab, col = col.axis),rot=90)
popViewport()

vpxl <- viewport(0.5,0.98,width=0.9,height = 0.05)
pushViewport(vpxl)
grid.text(label=x.name,0.5,0.5,gp=gpar(fontsize = 14*cex.lab,col=col.axis))
popViewport()


return(invisible(TRUE))
	
}



