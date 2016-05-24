GraphDistatisPartial <-
function(FS,PartialFS,axis1=1,axis2=2,constraints=NULL,item.colors=NULL,participant.colors=NULL,ZeTitle='Distatis-Partial',Ctr=NULL,color.by.observations=TRUE,nude=FALSE,lines=TRUE){
   
   
	if(is.null(participant.colors)){
	   	part.design <- diag(dim(PartialFS)[3])
		participant.colors <- as.matrix(createColorVectorsByDesign(part.design)$oc)
	}	
	if(is.null(item.colors)){
		item.design <- diag(dim(FS)[1])
		item.colors <- as.matrix(createColorVectorsByDesign(item.design)$oc)
	}	
	if(is.null(Ctr)){
		F2 <- FS**2
		SF2 <- apply(F2,2,sum)
		Ctr <- t( t(F2) / SF2)
	}	
	if(is.null(constraints)){
	    real.minimum <- min(c(FS,PartialFS))
    	real.maximum <- max(c(FS,PartialFS))
	    real.value <- max(c(abs(real.minimum),abs(real.maximum)))
   		#set up a constraints list
		constraints <- list(minx=-real.value, maxx=real.value,miny=-real.value,maxy=real.value)		
	}
	
	if(!color.by.observations){
		title.pass <- paste(ZeTitle,"Colored By Participants",sep=" ")	
		Compromise.Out <- GraphDistatisCompromise(FS,axis1=axis1,axis2=axis2,constraints=constraints,item.colors=item.colors,ZeTitle=title.pass,nude=nude,Ctr=Ctr)
		if(lines){
			for(i in 1:dim(PartialFS)[1]){
				to.plot <- t(PartialFS[i,,])		
				center.point <- FS[i,c(axis1,axis2)]
				center.rep <- matrix(center.point,dim(PartialFS)[3],2,byrow=TRUE)
				bound.mat <- rbind(center.rep,to.plot[,c(axis1,axis2)])
				bound.mat <- bound.mat[ as.vector(t(matrix(seq(1,nrow(bound.mat)),ncol=2))), ]
				points(bound.mat,type="l",lty=2,lwd=2,col="grey80")
			}
		}		
		for(i in 1:dim(PartialFS)[3]){
			to.plot <- PartialFS[,,i]
			rownames(to.plot) <- rep(unlist(dimnames(PartialFS)[3])[i],dim(PartialFS)[1])
			prettyPlot(to.plot,col=participant.colors[i,],dev.new=FALSE,axes=FALSE,new.plot=FALSE,x_axis=axis1,y_axis=axis2,display_names=!nude)
		}				
	}else{
		title.pass <- paste(ZeTitle,"Colored By Observations",sep=" ")
		Compromise.Out <- GraphDistatisCompromise(FS,axis1=axis1,axis2=axis2,constraints=constraints,item.colors=item.colors,ZeTitle=title.pass,nude=nude,Ctr=Ctr)
		if(lines){
			for(i in 1:dim(PartialFS)[1]){
				to.plot <- t(PartialFS[i,,])		
				center.point <- FS[i,c(axis1,axis2)]
				center.rep <- matrix(center.point,dim(PartialFS)[3],2,byrow=TRUE)
				bound.mat <- rbind(center.rep,to.plot[,c(axis1,axis2)])
				bound.mat <- bound.mat[ as.vector(t(matrix(seq(1,nrow(bound.mat)),ncol=2))), ]
				points(bound.mat,type="l",lty=2,lwd=2,col="grey80")
			}
		}		
		for(i in 1:dim(PartialFS)[1]){
			to.plot <- t(PartialFS[i,,])			
			prettyPlot(to.plot,col=item.colors[i,],dev.new=FALSE,axes=FALSE,new.plot=FALSE,x_axis=axis1,y_axis=axis2,display_names=!nude)
		}				
	}
	
	retour <- Compromise.Out
	retour$participant.colors=participant.colors
	return(retour)
}
