
circutile <- function(x, col = "grey", col.opt = list(),persp = FALSE){
	nd <- length(dim(x))
	stopifnot(nd == 2)
	
	if("line.col" %in% names(col.opt)){
		line.col <- col.opt$line.col
	}else{
		line.col <- alpha(1,0.2)
	}
	if("circle.col" %in% names(col.opt)){
		circle.col <- col.opt$circle.col
	}else{
		circle.col <- line.col
	}
	
	n <- nrow(x)
	m <- ncol(x)
	colv <- matrix(col,n,m, byrow=TRUE)
	
	grid.newpage()
	v1 <- viewport(0.5,0.5,width=0.9,height=0.9)
	pushViewport(v1)
	sc <- 0.8
	r <- rep((1:m)/m, each = n)*sc + 1-sc
	ang <- rep( (1:n)/n*2*pi, m)
	sapply(((1:m)/m)*sc + 1 - sc,function(z){
		grid.circle(0.5,0.5,r=z/2,gp=gpar(col=circle.col))
		
	})
	sapply((1:n)/n*2*pi,function(z){
				grid.lines( x=c( (1+sin(z)*(1/m*sc +1 - sc))/2,(1+sin(z))/2), y=c( (1+cos(z)*(1/m*sc +1 - sc))/2, (1+cos(z))/2),gp=gpar(col=line.col))
	})
	if(persp){
		z <- x/max(x)/(max(n,m)+1)/4*rep(seq(1-sc,1,sc/(m-1))/(1-sc),each=n)
	}else{
		z <- x/max(x)/(max(n,m)+1)/4
	}
	
	xc <-  (1+r*sin(ang))/2
	yc <-  (1+r*cos(ang))/2
	
	cc <- colv
	grid.circle(x = xc, y = yc, r = z,gp = gpar(fill=colv))
	#mapply( function(x,y,z,cc){
	#	grid.circle(x = (1+y*sin(x))/2, y = (1+y*cos(x))/2, r = z,gp = gpar(fill=cc))
	#}, x = ang, y = r, z = x/max(x)/4/(max(n,m)+1),cc=colv)
			
			
	return(invisible(TRUE))
}