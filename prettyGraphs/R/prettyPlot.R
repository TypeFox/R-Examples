prettyPlot <-
function(data_matrix,x_axis=1,y_axis=2,col=NULL,pch=NULL,cex=NULL,text.cex=NULL,pos=3,xlab="",ylab="",main="",display_names=TRUE,display_points=TRUE,constraints=NULL,contributionCircles=FALSE,contributions=NULL,axes=TRUE,fg.line.width=3,fg.type="l",fg.col="black",bg.line.width=1.5,bg.lty=3,bg.col = "black",flip=FALSE,asp=1,findBounds=TRUE,dev.new=TRUE,new.plot=TRUE){
	
	#I want to always send back colors and constraints.
	#I need a different type of checker here...
	if(is.null(col)){
		col <- colorVectorIsNull(data_matrix)$oc
	}
	col <- as.matrix(col)
	if(length(col)==1){
		col <- as.matrix(rep(col,nrow(data_matrix)))
	}else if(nrow(col)!=nrow(data_matrix)){
		col <- colorVectorIsNull(data_matrix)$oc
	}
	
	if(is.null(pch)){
		pch <- as.matrix(rep(21,nrow(data_matrix)))
	}
	pch <- as.matrix(pch)
	if(length(pch)==1){
		pch <- as.matrix(rep(pch,nrow(data_matrix)))
	}else if(nrow(pch)!=nrow(data_matrix)){
		pch <- as.matrix(rep(21,nrow(data_matrix)))
	}
	
	if(is.null(cex)){
		cex <- as.matrix(rep(1,nrow(data_matrix)))
	}
	cex <- as.matrix(cex)
	if(length(cex)==1){
		cex <- as.matrix(rep(cex,nrow(data_matrix)))
	}else if(nrow(cex)!=nrow(data_matrix)){
		cex <- as.matrix(rep(1,nrow(data_matrix)))
	}
	
	if(is.null(text.cex)){
		text.cex <- as.matrix(rep(0.8,nrow(data_matrix)))
	}
	
	text.cex <- as.matrix(text.cex)	
	if(length(text.cex)==1){
		text.cex <- as.matrix(rep(text.cex,nrow(data_matrix)))
	}else if(nrow(text.cex)!=nrow(data_matrix)){
		text.cex <- as.matrix(rep(0.8,nrow(data_matrix)))
	}
	
	#I only need constraints if I am making a new window.
	check.constraints <- minmaxHelper(data_matrix,data_matrix,x_axis,y_axis,findBounds=findBounds)	
	if(!is.null(constraints)){
	#if it is not null
		if(("minx" %in% names(constraints)) && ("maxx" %in% names(constraints)) && ("miny" %in% names(constraints)) && ("maxy" %in% names(constraints))){
			#and if it meets criteria, use it. This way, if FALSE, it should also go here.
			check.constraints <- list(minx=constraints$minx,miny=constraints$miny,maxx=constraints$maxx,maxy=constraints$maxy)
		}
	}
	constraints <- check.constraints	
	
	
	if(	(display_names==FALSE && display_points==FALSE) ){
		#For your health!
		print("Sorry, but you cannot have display_points and display_names set to FALSE. Nothing would have been plotted!")
		display_points<-TRUE
	}
	
	#this assumes you already have a device ready to go.
	if(dev.new){
		dev.new()	
	}
	
	if(new.plot){
		plot(c(0,0),c(0,0),type="n",col="white",axes=FALSE,xlab=xlab,ylab=ylab,ylim=c(constraints$miny,constraints$maxy),xlim=c(constraints$minx,constraints$maxx),main=main,asp=asp)
	}
			
	#make a new plot on a device.
	if(axes){		
		#determine axis points
		axis_list <- determineAxesPosition(constraints)		
		#plot the axes
		makeAxes(axis_list,fg.line.width= fg.line.width,fg.type= fg.type,fg.col= fg.col,bg.line.width= bg.line.width,bg.lty=bg.lty,bg.col = bg.col)
	}	
	#am I displaying points?
	if(display_points){
		plotPoints(data_matrix,col,x_axis,y_axis,cex=cex,pch=pch,contributionCircles=contributionCircles,contributions=contributions,flip=flip)
	}
	#am I displaying names?
	if(display_names){
		plotText(data_matrix,col,x_axis,y_axis,pos=pos,text.cex=text.cex,contributionCircles=contributionCircles,contributions=contributions)
	}		
	
	return(list(col=col,pch=pch,constraints=constraints))
}
