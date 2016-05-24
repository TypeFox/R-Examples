plotPoints <-
function(data_matrix,col,x_axis=1,y_axis=2,cex=1,pch=21,contributionCircles=FALSE,contributions=NULL,flip=FALSE){
		
	##the pch checks may not be necessary or may need to move up to prettyPlot.
	# if(is.null(pch)){
		# pch <- as.matrix(rep(21,nrow(data_matrix)))
	# }else{
		# pch<-as.matrix(pch)
	# }
	
	# if(nrow(pch)!=nrow(data_matrix)){
		# pch <- as.matrix(rep(21,nrow(data_matrix)))
		# print("pch dimension mismatch. Default pch used.")
	# }
	
	##pch, col, and cex better be good by the time they get here. the checking needs to occur in epGraphs & prettyPlot.
	
	pch.with.bg <- c(21:25)
	pch.with.col <- c(1:20)	
	bg.indices <- which(pch %in% pch.with.bg)
	col.indices <- which(pch %in% pch.with.col)
	
	#cex computation.
	if(contributionCircles){
		if( !is.null(contributions) && !(nrow(as.matrix(contributions))!=nrow(data_matrix)) ){
			thesecontributions <- rowSums(contributions[,c(x_axis,y_axis)]^1.5)
			thesecontributions <- (thesecontributions / ((max(thesecontributions) - min(thesecontributions)) / (2 - 0))) + 1			
			thesecontributions[which(thesecontributions > 3)] <- 3
			cex <- as.matrix(cex * thesecontributions)
		}else{
			print("Contributions dimension mismatch. Defaulting to cex only.")
		}
	}
	
	##plotting happens here.
	if(length(bg.indices) > 0){
		if(flip){
			points(data_matrix[bg.indices,x_axis],data_matrix[bg.indices,y_axis],col=col[bg.indices,],pch=pch[bg.indices,],bg="black",cex=cex[bg.indices,])			
		}else{
			points(data_matrix[bg.indices,x_axis],data_matrix[bg.indices,y_axis],col="black",pch=pch[bg.indices,],bg=col[bg.indices,],cex=cex[bg.indices,])
		}
	}
	
	if(length(col.indices) > 0){
		points(data_matrix[col.indices,x_axis],data_matrix[col.indices,y_axis],col=col[col.indices,],pch=pch[col.indices,],cex=cex[col.indices,])		
	}
	
	
}
