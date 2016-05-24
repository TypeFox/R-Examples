makeAxes <-
function(axis_position_list,fg.line.width=3,fg.type="l",fg.col="black",bg.line.width=1.5,bg.lty=3,bg.col = "black"){
	#This should make a blank plot window with custom placed axes and labels.
	#dev.new()
	#not sure yet...
#	plot(c(0,0),c(0,0),type="n",col="white",axes=FALSE,xlab=xlab,ylab=ylab,ylim=c(min_max_list$miny,min_max_list$maxy),xlim=c(min_max_list$minx,min_max_list$maxx),main=main,asp=asp)
    abline(v = 0, lty = bg.lty, lwd= bg.line.width, col = bg.col)
    abline(h = 0, lty = bg.lty, lwd= bg.line.width, col = bg.col)
    
	points(axis_position_list$xx,axis_position_list$xy,type=fg.type,col= fg.col,lwd=fg.line.width)
	points(axis_position_list$yx,axis_position_list$yy,type=fg.type,col= fg.col,lwd=fg.line.width)
}
