plot.rrs <- function(x, ..., display = c("heat","contour","persp","3d"), show.WIN = TRUE, tolerance.matrix = NULL, tol.opt = list(raise = 0.01, col = "black", levels = 0.05, lty = 1, lwd = 1)){
	extras <- list(...)
	if(is.null(extras$xlab)) extras$xlab <- ""
	if(is.null(extras$ylab)) extras$ylab <- ""
	if(is.null(extras$zlab)&&(display=="persp"||display=="3d")) extras$zlab <- ""
	if(is.null(extras$main)) extras$main <- ""
	if(length(display)==4||display=="heat"){
		image <- im(t(x$rsM),xcol=x$f$X,yrow=x$f$Y)[x$f$WIN,drop=F]
		if(is.null(list(...)$axes)) plot(image,...,axes=T)
		else plot(image,...)
		if(show.WIN) plot(x$f$WIN,add=T,lwd=2)
		if(!is.null(tolerance.matrix)){
			gr <- nrow(tolerance.matrix)
			
			if(is.null(tol.opt$col)) tol.opt$col <- "black"
			if(is.null(tol.opt$levels)) tol.opt$levels <- 0.05
			if(is.null(tol.opt$lty)) tol.opt$lty <- 1
			if(is.null(tol.opt$lwd)) tol.opt$lwd <- 1
			
			contour(seq(x$f$WIN$xrange[1],x$f$WIN$xrange[2],length=gr),seq(x$f$WIN$yrange[1],x$f$WIN$yrange[2],length=gr),tolerance.matrix,add=T,drawlabels=F,levels=tol.opt$levels,lty=tol.opt$lty,lwd=tol.opt$lwd,col=tol.opt$col)
		}
	} else if(display=="contour"){
		contour(x$f$X,x$f$Y,x$rsM,...)
		if(show.WIN) plot(x$f$WIN,add=T,lwd=2)
	} else if(display=="persp"){
		if(!is.null(extras$col)){
			if(length(extras$col)>1){
				zfacet <- x$rsM[-1,-1]+x$rsM[-1,-ncol(x$rsM)]+x$rsM[-nrow(x$rsM),-1]+x$rsM[-nrow(x$rsM),-ncol(x$rsM)]
				facetcol <- cut(zfacet,length(extras$col))
				extras$col <- extras$col[facetcol]
			}
		}
		do.call("persp",c(list(x=x$f$X,y=x$f$Y,z=x$rsM), extras))
		
	} else if(display=="3d"){
		if(!is.null(extras$col)){
			if(length(extras$col)>1){
				extras$col <- assignColors(extras$col,as.matrix(im(t(x$rsM),xcol=x$f$X,yrow=x$f$Y)))
			}
		}
		do.call("persp3d", c(list(x=x$f$X,y=x$f$Y,z=x$rsM), extras)) 
		if(show.WIN){
			w <- x$f$WIN
			wv <- vertices(w)
			tempdf <- as.data.frame(wv)
			gridLocs <- rep(NA,nrow(tempdf))
						
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
		if(!is.null(tolerance.matrix)){
			
			gr <- nrow(tolerance.matrix)
			xt <- sort(rep(seq(x$f$WIN$xrange[1],x$f$WIN$xrange[2],length=gr),gr))
			yt <- rep(seq(x$f$WIN$yrange[1],x$f$WIN$yrange[2],length=gr),gr)
			tempmat <- matrix(c(xt,yt),gr*gr,2)
			
			tempmask <- as.im(x)
			grdcoo <- nearest.raster.point(x=tempmat[,1],y=tempmat[,2],w=tempmask)
			
			rhot <- rep(NA,nrow(tempmat))
			for(i in 1:nrow(tempmat)) rhot[i] <- tempmask$v[grdcoo$row[i],grdcoo$col[i]]
			
			if(is.null(tol.opt$raise)) tol.opt$raise <- 0.01
			if(is.null(tol.opt$col)) tol.opt$col <- "black"
			if(is.null(tol.opt$levels)) tol.opt$levels <- 0.05
			if(is.null(tol.opt$lwd)) tol.opt$lwd <- 1
			
			tol3d(seq(x$f$WIN$xrange[1],x$f$WIN$xrange[2],length=gr),seq(x$f$WIN$yrange[1],x$f$WIN$yrange[2],length=gr),tolerance.matrix,rhot,tol.opt$levels,tol.opt$raise,tol.opt$col,tol.opt$lwd)
		}
	}
}
