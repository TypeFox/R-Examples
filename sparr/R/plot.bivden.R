plot.bivden <- function(x, ..., display = c("heat","contour","persp","3d"), show.WIN = TRUE){
	extras <- list(...)
	if(is.null(extras$xlab)) extras$xlab <- ""
	if(is.null(extras$ylab)) extras$ylab <- ""
	if(is.null(extras$zlab)&&(display=="persp"||display=="3d")) extras$zlab <- ""
	if(is.null(extras$main)) extras$main <- ""
	if(length(display)==4||display=="heat"){
		image <- im(t(x$Zm),xcol=x$X,yrow=x$Y)[x$WIN,drop=F]
		if(is.null(list(...)$axes)) plot(image,...,axes=T)
		else plot(image,...)
		if(show.WIN) plot(x$WIN,add=T,lwd=2)
	} else if(display=="contour"){
		contour(x$X,x$Y,x$Zm,...)
		if(show.WIN) plot(x$WIN,add=T,lwd=2)
	} else if(display=="persp"){
		if(!is.null(extras$col)){
			if(length(extras$col)>1){
				zfacet <- x$Zm[-1,-1]+x$Zm[-1,-ncol(x$Zm)]+x$Zm[-nrow(x$Zm),-1]+x$Zm[-nrow(x$Zm),-ncol(x$Zm)]
				facetcol <- cut(zfacet,length(extras$col))
				extras$col <- extras$col[facetcol]
			}
		}
		do.call("persp",c(list(x=x$X,y=x$Y,z=x$Zm), extras))

	} else if(display=="3d"){
		if(!is.null(extras$col)){
			if(length(extras$col)>1){
				extras$col <- assignColors(extras$col,as.matrix(im(t(x$Zm),xcol=x$X,yrow=x$Y)))
			}
		}
		do.call("persp3d", c(list(x=x$X,y=x$Y,z=x$Zm), extras)) 
		if(show.WIN){
			w <- x$WIN
			wv <- vertices(w)
			tempdf <- as.data.frame(wv)
			tempmask <- as.im(x)
			grdcoo <- nearest.raster.point(x=tempdf[,1],y=tempdf[,2],w=tempmask)
			
			zgrd <- rep(NA,nrow(tempdf))
			for(i in 1:nrow(tempdf)) zgrd[i] <- tempmask$v[grdcoo$row[i],grdcoo$col[i]]
			
			zgrdw <- c(zgrd,zgrd[1])
			xl <- c(wv$x,wv$x[1])[-which(is.na(zgrdw))]
			yl <- c(wv$y,wv$y[1])[-which(is.na(zgrdw))]
			zgrdw <- zgrdw[-which(is.na(zgrdw))]

			lines3d(xl,yl,zgrdw,lwd=4)
		}
	}
}
