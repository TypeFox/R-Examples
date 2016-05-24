correlationPlotter <-
function(data_matrix,factor_scores,x_axis=1,y_axis=2,col=NULL,pch=NULL,xlab=NULL,ylab=NULL,main="",asp=1,dev.new=TRUE){
	if(nrow(data_matrix)==nrow(factor_scores)){
		loadings <- cor(data_matrix,factor_scores)
	}
	else if(ncol(data_matrix)==nrow(factor_scores)){
		loadings <- cor(t(data_matrix),factor_scores)
	}else{
		print("Dimension mismatch. Please check data_matrix and factor_scores.")
		return(NULL)	
	}
	loadings <- replace(loadings,is.na(loadings),0)

	if(dev.new){
		dev.new()
	}
	
	if(!is.null(xlab) && !is.null(ylab)){
		plotCircle(xlab=xlab,ylab=ylab,main=main,asp=asp)	
	}else{
		plotCircle(xlab=paste("Component ",x_axis,sep=""),ylab=paste("Component ",y_axis,sep=""),main=main,asp=asp)
	}
	
	if(is.null(col)){
		col <- colorVectorIsNull(loadings)$oc
	}
	if(is.null(pch)){
		pch <- as.matrix(rep(21,nrow(loadings)))
	}	


	os <- cbind(rep(0,nrow(loadings)+1),rep(0,nrow(loadings)+1))
	new.mat <- matrix(0,(nrow(loadings)*2)+1,2)
	new.mat[seq(1,nrow(new.mat),2),] <- os
	new.mat[seq(2,nrow(new.mat),2),] <- loadings[,c(x_axis,y_axis)]
	points(new.mat,type="l",col="black")
	
	prettyPlot(loadings,col=col,display_names=TRUE,display_points=TRUE,pch=pch,x_axis=x_axis,y_axis=y_axis,axes=FALSE,dev.new=FALSE,new.plot=FALSE)
}
