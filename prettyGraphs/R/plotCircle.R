plotCircle <-
function(xlab="",ylab="",main="",axis_line_width=3,asp=1){

	j = seq(pi/128,2*pi,by=pi/128)
	coords = cbind(cos(j),sin(j))
	coords = rbind(coords,coords[1,])
	#an override because I could not trace the real problem...
	constraints <- list(minx=-1.1,miny=-1.1,maxx=1.1,maxy=1.1)
	axis_list <- determineAxesPosition(constraints)
	#
	#makeNewPlotWindow(axis_list,constraints,xlab,ylab,main,axis_line_width,asp=asp)
	plot(c(0,0),c(0,0),type="n",col="white",axes=FALSE,xlab=xlab,ylab=ylab,ylim=c(constraints$miny,constraints$maxy),xlim=c(constraints$minx,constraints$maxx),main=main,asp=asp)	
	makeAxes(axis_list,axis_line_width)	
	points(coords,col="black",type='l')
	points(0,0,col="black",pch=20)
}
