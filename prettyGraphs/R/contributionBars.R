contributionBars <-
function(factor_scores,contributions,x_axis=1,y_axis=2,col=NULL,main=NULL,upper='steelblue4',lower='firebrick2',threshold=0,sortContributions=TRUE,pretty=FALSE,show.bg.bars=FALSE){

	##run the checks.
	factor_scores <- as.matrix(factor_scores)
	contributions <- as.matrix(contributions)
	if(sum(rownames(factor_scores)==rownames(contributions)) != nrow(factor_scores)){
		rownames(contributions) <- rownames(factor_scores)
	}
	#multiply the sign of the data by the contributions
		#stupid rounding somewhere...
#	print(!(sum(colSums(contributions)>0.999)==ncol(contributions)))
#	print(colSums(contributions)>0.999 | colSums(contributions)==0)
#	print(sum(colSums(contributions)>0.999 | colSums(contributions)==0)!=ncol(contributions))
#	pause()
#	if(!(sum(colSums(contributions)>0.999)==ncol(contributions))){
	if(sum(colSums(contributions)>0.999 | colSums(contributions)==0)!=ncol(contributions)){
		stop(paste("Contributions do not sum to 1",which( (colSums(contributions)==1)==FALSE),sep=" "))
	}
	contributions_with_signs <- sign(factor_scores) * contributions
	
	if(threshold == 0){
		threshold <- 1/dim(contributions)[1]
	} else if(threshold >= 1){
		threshold <- 1/dim(contributions)[1]
	}
	
	resetCols <- FALSE
	if(is.null(col)){
		resetCols <- TRUE
	}

	#contribute those bars.
	dev.new()
	par(mfrow=c(1,2))
	flipFlag<-TRUE
	axes <- c(x_axis,y_axis)
	
	for(i in 1:length(axes)){
		if(resetCols){
			col <- matrix("gray",nrow(contributions_with_signs),1)
			col[which(contributions_with_signs[,i] <= -threshold),1] <- "darkseagreen"
			col[which(contributions_with_signs[,i] >= threshold),1] <- "plum4"
		}		
		if(sortContributions){
			ordered_inds <- order(contributions_with_signs[,axes[i]],decreasing=FALSE)
		}else{
			ordered_inds <- 1:dim(contributions_with_signs)[1]
		}
		ordered <- contributions_with_signs[ordered_inds,axes[i]]	
		ordered_colors <- col[ordered_inds]
		
		#draw a line across the thresholds.
		if(!flipFlag){
			#horizontal cut lines
			if(!pretty){
				barplot(ordered,col=ordered_colors,ylim=c(-1.1,1.1),axes=TRUE,horiz=flipFlag,sub=paste("Component ",axes[i],sep=""))
				abline(h=0,col="black")				
			}else{
				prettyBars(ordered,fg.col=ordered_colors,axis.lims=c(-1.1,1.1),horiz=flipFlag,dev.new=FALSE,bg.lims=c(-1,1),show.bg.bars=show.bg.bars)
			}
			abline(h=threshold,col=upper,lty=2)
			abline(h=-threshold,col=lower,lty=2)
		}else{
			#vertical cut lines		
			if(!pretty){	
				barplot(ordered,col=ordered_colors,xlim=c(-1.1,1.1),axes=TRUE,horiz=flipFlag,sub=paste("Component ",axes[i],sep=""))
				abline(v=0,col="black")								
			}else{
				prettyBars(ordered,fg.col=ordered_colors,axis.lims=c(-1.1,1.1),horiz=flipFlag,dev.new=FALSE,bg.lims=c(-1,1),show.bg.bars=show.bg.bars)				
			}
			abline(v=threshold,col=upper,lty=2)
			abline(v=-threshold,col=lower,lty=2)			
		}
		flipFlag<-!flipFlag		
		#make a pie chart of contributions
		#later!
	}

	if(is.null(main)){		
		mtext("Contributions to variance",side=3,outer=TRUE,line=-3)
	}else{
		mtext(main,side=3,outer=TRUE,line=-3)
	}

}
