#wonderful stuff from here: http://martin.ankerl.com/2009/12/09/how-to-create-random-colors-programmatically/


prettyGraphsHSVColorSelection <- function(n.colors=1,offset=NULL,h=13,s=0.75,v=0.75){

	if(is.null(offset)){
		offset <- abs((1-sqrt(5))/2)
	}
	these.cols <- c()
	for(i in 1:n.colors){
		h <- h + offset
		h <- h %% 1
		these.cols <- c(these.cols,hsv(h,s,v))
	}
	return(as.matrix(these.cols))
}

#plot(1:ncolors,1:ncolors,col=these.cols,cex=3,pch=20)