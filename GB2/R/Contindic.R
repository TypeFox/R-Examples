contindic.gb2 <- function(resol, shape1, shape21, shape22, shape31, shape32, fn, title, table = FALSE){

	pp <- round(seq(shape21,shape22,length.out=resol),digits=2)
	qq <- round(seq(shape31,shape32,length.out=resol),digits=2)
	d <- pp %o% qq
	for (i in 1:resol){
		for (j in 1:resol){
		d[i,j] <- fn(shape1,pp[i],qq[j])
		}
	}
  dlim = range(d, finite = TRUE)
	contour(pp, qq, d, levels = pretty(seq(dlim[1], dlim[2], length.out=17)),
	xlab = "p", ylab = "q", main = paste("a =", as.character(shape1)), cex=1.8)
	box()
  mtext(title,line=0.5)
	if(table){
		d <- rbind(p=pp, d) 
		d <- round(cbind(c(NA,qq), d), digits=2)	
	print(paste("a =", as.character(shape1)), quote=FALSE)
	print(d)
	}
}
