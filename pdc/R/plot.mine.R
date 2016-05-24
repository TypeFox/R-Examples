plot.mine <- function(x, normalize=FALSE,type="image", mark.optimum=TRUE, col=heat.colors(12), ...) {

x.data <- x$t.range
y.data <- x$m.range

if ((length(x.data) == 1) | (length(y.data) == 1)) {
	
 if (length(x.data)==1) {
 	xlabel <- "embedding dimension"
 	rng <- y.data
 	r <- x$m
 	idx <- 2
 } else {
 	idx <- 1
 	rng <- x.data
 	r <- x$t
 	xlabel <- "time delay"
 }
	
 pch <- rep(22, length(rng))
 
 if (mark.optimum) {
  pch[ which(rng==r)] <- 15
 } 

 plot(x$entropy.values[,idx], x$entropy.values[,3], type="b", pch=pch,
 	xlab=xlabel, ylab="entropy",...)
 if (mark.optimum) {
 	abline(v = r, lty=3)
 }

	
	
} else {

 z <- matrix(x$entropy.values[,3],nrow=length(x.data),byrow=T)
 if (normalize)
	z <- (z-min(z))/(max(z)-min(z))

# if (type == "interp.image") {
#	 z <- interp(x$entropy.values[,1],x$entropy.values[,2],x$entropy.values[,3])
#	 image(z, xlab="time delay",ylab="embedding dimension",main="Entropy Heuristic")
# } else 
 if (type=="contour") {
	contour(x.data,y.data,z, ...)
 #} else if (type=="image") {
  } else {
	image(x.data,y.data,z, xlab="time delay",ylab="embedding dimension",main="Entropy Heuristic" , col=col, ...)
 }
 
 
 # add image scale
#	if (scale) {
		
	#}
	
 if (mark.optimum) {
 	points( x$t, x$m, pch=4, cex=3)
 }

}

invisible()

}