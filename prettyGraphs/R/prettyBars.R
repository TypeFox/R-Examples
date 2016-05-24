prettyBars <- function(data,axis=1,cex.names=0.5,fg.col=NULL,axis.lims=NULL,show.bg.bars=FALSE,threshold.line=FALSE,main="",bg.border="white",bg.col=NULL,bg.lims=NULL,sort.data=TRUE,horiz=TRUE,dev.new=TRUE){
	
	#I want data to be a matrix, for now.
	data <- as.matrix(data)
	
	#small fix... 
	max.value <- max(abs(data))
	if(is.null(axis.lims)){
		axis.lims <- c(-max.value * 1.2,max.value * 1.2)		
		if(!is.null(bg.lims) && length(bg.lims)==2){
			if(axis.lims[1] > bg.lims[1]){
				axis.lims[1] <- bg.lims[1]
			}
			if(axis.lims[2] < bg.lims[2]){
				axis.lims[2] <- bg.lims[2]
			}
		}			
	}
	
	sign.values <- sign(data[,axis])
	poss <- sum(which(sign.values>= 0))
	negs <- sum(which(sign.values==-1))
		
	if(is.null(fg.col)){
		#create colors
		fg.col <- matrix("plum4",nrow(data),1)
		fg.col[which(sign.values==-1),1] <- "darkseagreen"		
	}
	
	if(sort.data){
		order.inds <- order(data[,axis])
		data <- data[order.inds,]
		fg.col <- as.matrix(fg.col[order.inds,])
		##there should be a faster way...	
		sign.values <- sign(data[,axis])
		poss <- sum(which(sign.values>= 0))
		negs <- sum(which(sign.values==-1))		
	}
	
	if(dev.new){
		dev.new()
	}
	if(is.null(bg.col)){
		bg.col <- c("gray","gray")
	}	
	if(show.bg.bars){
		sign.values.mat <- rbind(sign.values, sign.values)
		#update this to allow for a sum all inputted values. or at least some update of it
		if(is.null(bg.lims) || length(bg.lims) != 2){
			sals.bound <- rbind(abs(data[,axis]), max.value-abs(data[,axis]))
			
		}else{
			sals.bound <- rbind(rep(bg.lims[2],nrow(data)),rep(bg.lims[1],nrow(data)))
		}
		barplot.values <- sign.values.mat * sals.bound		
		rownames(barplot.values) <- NULL	
		colnames(barplot.values) <- NULL
		if(horiz){		
			bp.bg <- barplot(barplot.values,beside=FALSE,horiz=horiz,xlim=axis.lims,axes=FALSE,border=bg.border,col=bg.col,main=main)
		}else{
			bp.bg <- barplot(barplot.values,beside=FALSE,horiz=horiz,ylim=axis.lims,axes=FALSE,border=bg.border,col=bg.col,main=main)
		}
	}

	data.copy <- data
	rownames(data.copy) <- NULL
	colnames(data.copy) <- NULL
	if(horiz){
		bp.cols <- barplot(data.copy[,axis],col=fg.col,horiz=horiz,xlim=axis.lims,axes=FALSE,border=fg.col,add=show.bg.bars,main=main)
		abline(v=0,lty=5,lwd=2)
		if(threshold.line && (!is.null(bg.lims) && length(bg.lims)==2)){
			abline(v=bg.lims,lwd=2,col=bg.col,lty=4)
		}
		if(poss){
			text(data[which(sign.values >= 0),axis], bp.cols[which(sign.values >= 0)], rownames(data)[which(sign.values >= 0)],cex=cex.names,adj=-0.1,col=fg.col[which(sign.values >= 0),])
		}
		if(negs){
			text(data[which(sign.values ==-1),axis], bp.cols[which(sign.values ==-1)], rownames(data)[which(sign.values ==-1)],cex=cex.names,adj=1.1,col=fg.col[which(sign.values == -1),])
		}
	}else{
		bp.cols <- barplot(data.copy[,axis],col=fg.col,horiz=horiz,ylim=axis.lims,axes=FALSE,border=fg.col,add=show.bg.bars,main=main)
		abline(h=0,lty=5,lwd=2)
		if(threshold.line && (!is.null(bg.lims) && length(bg.lims)==2)){
			abline(h=bg.lims,lwd=2,col=bg.col,lty=4)
		}		
		if(poss){
			text(bp.cols[which(sign.values >= 0)], data[which(sign.values >= 0),axis],  rownames(data)[which(sign.values >= 0)],cex=cex.names,adj=-0.1,srt=90,col=fg.col[which(sign.values >= 0),])
		}
		if(negs){
			text(bp.cols[which(sign.values ==-1)], data[which(sign.values ==-1),axis], rownames(data)[which(sign.values ==-1)],cex=cex.names,adj=1.1,srt=90,col=fg.col[which(sign.values == -1),])		
		}
	}
	
	return(bp.cols)
}