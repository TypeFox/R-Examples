GraphDistatisAll <-
function(FS,PartialFS,FBoot,RvFS,axis1=1,axis2=2,constraints=NULL,item.colors=NULL,participant.colors=NULL,ZeTitleBase=NULL,nude=FALSE,Ctr=NULL,RvCtr=NULL,color.by.observations=TRUE,lines=TRUE,lwd=3.5,ellipses=TRUE,fill=TRUE, fill.alpha = .27, percentage=0.95){

	if(is.null(participant.colors)){
	   	part.design <- diag(dim(PartialFS)[3])
		participant.colors <- as.matrix(createColorVectorsByDesign(part.design)$oc)
	}	
	if(is.null(item.colors)){
		item.design <- diag(dim(FS)[1])
		item.colors <- as.matrix(createColorVectorsByDesign(item.design)$oc)
	}
	if(is.null(ZeTitleBase)){
		ZeTitleBase <- 'Distatis'
	}
	if(is.null(Ctr)){
		F2 <- FS**2
		SF2 <- apply(F2,2,sum)
		Ctr <- t( t(F2) / SF2)
	}
	if(is.null(RvCtr)){
		RvF2 <- RvFS**2
		RvSF2 <- apply(RvF2,2,sum)
		RvCtr <- t( t(RvF2) / RvSF2)
	}	
	if(is.null(constraints)){
    	#First, compute the constraints
	    real.minimum <- min(c(PartialFS,FS,FBoot))
    	real.maximum <- max(c(PartialFS,FS,FBoot))
	    real.value <- max(c(abs(real.minimum),abs(real.maximum)))
   		#set up a constraints list
		constraints <- list(minx=-real.value, maxx=real.value,miny=-real.value,maxy=real.value)
   }	
	
	compromise.out <- GraphDistatisCompromise(FS=FS,axis1=axis1,axis2=axis2,constraints=constraints,item.colors=item.colors,ZeTitle=paste(ZeTitleBase,"Compromise",sep=" "),nude=nude,Ctr=Ctr)

	partial.out <- GraphDistatisPartial(FS=FS,PartialFS=PartialFS,axis1=axis1,axis2=axis2,constraints=constraints,item.colors=item.colors,participant.colors=participant.colors,ZeTitle=paste(ZeTitleBase,"Partial",sep=" "), color.by.observations=color.by.observations,lines=lines)

	boot.out <- GraphDistatisBoot(FS=FS,FBoot=FBoot,axis1=axis1, axis2=axis2, item.colors=item.colors, ZeTitle=paste(ZeTitleBase,"Bootstrap",sep=" "), constraints=constraints, nude = nude, Ctr=Ctr, lwd = lwd, ellipses=ellipses, fill = fill, fill.alpha = fill.alpha, percentage=percentage)	

	rv.map.out <- GraphDistatisRv(RvFS,axis1=axis1,axis2=axis2,ZeTitle=paste(ZeTitleBase,"Rv Map",sep=" "), participant.colors = participant.colors, nude=nude,RvCtr=RvCtr)

	return(list(constraints=constraints,item.colors=item.colors,participant.colors=participant.colors))
	
}
