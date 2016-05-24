colormap <-
function(modgamobj, map=NULL, add=F, contours="none", mapmin=NULL, mapmax=NULL, 
	arrow=T, axes=F, ptsize=0.9, alpha=0.05) {
	if (!is.null(modgamobj$OR)) fitvals=modgamobj$OR else fitvals=modgamobj$fit 
	results = cbind(modgamobj$grid,fitvals)
	if (!is.null(modgamobj$pointwise)) results = cbind(results,modgamobj$pointwise)
	if (is.null(mapmin)) mapmin=min(results[,3])
	if (is.null(mapmax)) mapmax=max(results[,3])
	if (!is.null(map)) {
		if (add==T) {add=F; warning("add=T ignored because the map argument was used")}  
		if (class(map)=="map") maprange = map$range else
		if (class(map)=="SpatialPolygonsDataFrame") {
			mapcenter = apply(bbox(map),1,mean)
			center2edges = c(-1,1)*max(apply(bbox(map),1,diff))/2
			maprange = c(mapcenter[1]+center2edges,mapcenter[2]+center2edges)
		} else
		maprange = c(Inf,-Inf,Inf,-Inf)
	}
	dataXmin=min(results[,1],if(!is.null(map)) maprange[1])
	dataXmax=max(results[,1],if(!is.null(map)) maprange[2])
	dataYmin=min(results[,2],if(!is.null(map)) maprange[3])
	dataYmax=max(results[,2],if(!is.null(map)) maprange[4])
	offsetX=(dataXmax-dataXmin)/5
	offsetY=(dataYmax-dataYmin)/5
	qu = seq(mapmin,mapmax,length=2251)
	cp = rainbow(2252,start=0,end=0.66)
	col.seq = rev(cp)
	grad = cut(results[,3],breaks=c(0,qu,Inf),labels=F)
	par(xpd=TRUE)
	if (add==T) points(modgamobj$grid,col=col.seq[grad],pch=15,cex=ptsize) else
	if (axes==F) {
		par(mai=c(0,0.35,0.7,0.35))
		plot(modgamobj$grid,col=col.seq[grad],pch=15,cex=ptsize,type="p",
			 xlim=c(dataXmin-0.4*offsetX,dataXmax+0.4*offsetX),
			 ylim=c(dataYmin-0.8*offsetY,dataYmax),ann=F,axes=F)
	   }  else
	if (axes==T) {
#		par(mai=c(1.52,0.82,0.82,1.42))
		par(mai=c(1,0.82,1,0.5))
		plot(modgamobj$grid,col=col.seq[grad],pch=15,cex=ptsize,type="p",
			xlab="", xaxt="n")
		axis(3)
		mtext(names(modgamobj$grid[1]), side=3, line=3)
	}
	if (!is.null(map)) {
		if (class(map)=="SpatialPolygonsDataFrame") plot(map, border=gray(0.3), add=T) else {
			if (class(map)!="map") warning("map class not recognized--attempting generic plot")
			lines(map,col=gray(0.3),type="l")
		}
	}
	
	# CONTOUR LINES / ISOBOLES
	if (contours=="permrank" && dim(results)[2] < 4) {
		warning("permrank contours omitted because permute=0 in modgam")
		contours = "none"
	}
	if (contours!="none") {		# if contour lines should be drawn
		if (contours=="response") concol=3 else 
			if (contours=="permrank") concol=4  else
			stop("contours type not recognized")
		df=cbind(results[(order(results[,2])),1:2],results[(order(results[,2])),concol])
		cx=as.matrix(unique(df[(order(df[,1])),1]))
		cy=as.matrix(unique(df[,2]))
		cz=as.matrix(reshape(df,v.names=names(df)[3],idvar=names(df)[1],
			timevar=names(df)[2],direction="wide")[,-1])
		cz=cz[order(unique(df[,1])),]	
		if (contours=="response") contour(x=cx, y=cy, z=cz, add=TRUE)
		if (contours=="permrank") contour(x=cx, y=cy, z=cz, add=TRUE,
		 	levels=c(alpha/2,1-alpha/2), drawlabels=FALSE, lwd=3, lty="twodash")
	}
		
	# LEGEND
	if (!is.null(modgamobj$OR)) leglab = "Odds Ratios" else leglab = "Predicted Values"
	if (!is.null(modgamobj$m)) leglab = paste(modgamobj$m,leglab)
    fY = 0.5 + 0.5*(axes==T)
	ypos = dataYmin-fY*offsetY
	len = (dataXmax-dataXmin)*7/12
	points(dataXmin+(1:2252*len/2252),rep(ypos,2252),cex=1.6,col=cp[2253-1:2252],pch=15)
	if (!is.null(modgamobj$OR))  points(dataXmin+len*(1-mapmin)/(mapmax-mapmin),
		ypos,cex=.8,col="black", pch="|", lwd=2)  # only add line at 1 for ORs
	text(x=dataXmin,y=ypos,pos=1,labels=format(round(mapmin,2),digits=3),cex=.8)
	text(x=dataXmin+len,y=ypos,pos=1,labels=format(round(mapmax,2),digits=3),cex=.8)
	if (!is.null(modgamobj$OR)) text(x=dataXmin+len*(1-mapmin)/(mapmax-mapmin),
		y=ypos,pos=1,labels="1.00",cex=.8)	# only add line at 1 for ORs
	text(x=dataXmin+len/2,y=ypos,pos=3,labels=leglab,cex=1) 

	# NORTH ARROW
	if (arrow) {
		points(dataXmax,dataYmin-0.65*offsetY,cex=1.2,col="black", pch="|", lwd=2)
		points(dataXmax-0.01*offsetX,ypos,col="black", bg="black", pch=24)
		text(x=dataXmax,y=ypos,"N",cex=1,pos=3)
	}

	# SCALE BAR
	if (is.null(map)==F & class(map)=="SpatialPolygonsDataFrame") {
		d = 0.17*(dataXmax - dataXmin)
		d = signif(d,1)
		leftedge = dataXmin+1.125*len
		lines(c(leftedge,leftedge+d),c(ypos,ypos),col="black",lwd=3)
		points(leftedge,ypos,cex=1.2,col="black", pch="|", lwd=2)
		points(leftedge+d,ypos,cex=1.2,col="black",pch="|",lwd=2)
		if (d>=1000) {scdist=d/1000; units="km"} else {scdist=d; units="m"}
		sclab = paste(as.character(format(scdist,digits=2)),units)
		text(x=leftedge+0.5*d,y=ypos,sclab,cex=0.8,pos=1)
	}
}
