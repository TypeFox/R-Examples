#' @importFrom lattice lsegments
#' @importFrom lattice ltext
#' @importFrom rgl lines3d
#' @importFrom graphics segments
#' @importFrom graphics text
addarrow <- function(x0,y0,len,bearing,headlen=0.2*len,headangle=25,Nlabel=TRUE,
	addto=c('rgl','plot','lattice'),...){

	addto <- match.arg(addto)
	
	bearing <- bearing * pi/180
	headangle <- headangle * pi/180
	
	x1 <- x0 - len*cos(bearing)
	y1 <- y0 - len*sin(bearing)

	x2 <- x1 + headlen*cos(headangle + bearing)
	y2 <- y1 + headlen*sin(headangle + bearing)
	
	x3 <- x1 + headlen*cos(headangle - bearing)
	y3 <- y1 - headlen*sin(headangle - bearing)

	if(addto=="rgl"){
	  lines3d(x = c(x0,x1), y=c(y0,y1), z=c(0,0), ...)
	  lines3d(x = c(x1,x2), y=c(y1,y2), z=c(0,0), ...)
	  lines3d(x = c(x1,x3), y=c(y1,y3), z=c(0,0), ...)
	}
	if(addto=="lattice"){
	 lsegments(x0,y0,x1,y1,...)	
	 lsegments(x1,y1,x2,y2,...)
	 lsegments(x1,y1,x3,y3,...)
	}
	if(addto=="plot"){
	 segments(x0,y0,x1,y1,...)	
	 segments(x1,y1,x2,y2,...)
	 segments(x1,y1,x3,y3,...)
	}
	if(Nlabel){
		xN <- x0 + len*0.1*cos(bearing)
		yN <- y0 + len*0.1*sin(bearing)
		if(addto == "rgl")text3d(x=xN,y=yN,z=0,texts="N")
		if(addto == "plot")text(x=xN,y=yN,labels="N")
		if(addto == "lattice")ltext(x=xN,y=yN,labels="N")	
	}
}







