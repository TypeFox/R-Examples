GraphDistatisRv <-
function(RvFS,axis1=1,axis2=2,ZeTitle= 'Distatis-Rv Map', participant.colors = NULL, nude=FALSE,RvCtr=NULL){

	if(is.null(participant.colors)){	  
		participant.design <- diag(dim(RvFS)[1])
		participant.colors <- as.matrix(createColorVectorsByDesign(participant.design)$oc)
	}
	
	LeMinimum =  apply(RvFS,2,min)
	LeMaximum =   apply(RvFS,2,max)
	
	petitx =  min(c(0,LeMinimum[axis1]));grandx = LeMaximum[axis1]
	petity =  LeMinimum[axis2];grandy = LeMaximum[axis2]
	fudgeFact_H = (grandx-petitx)/9
	fudgeFact_V = (grandy-petity)/9
	constraints <- list(minx=petitx-fudgeFact_H,maxx=grandx+fudgeFact_H,miny=petity-fudgeFact_V,maxy=grandy+fudgeFact_V)
	  
	if(is.null(RvCtr)){
		RvF2 <- RvFS**2
		RvSF2 <- apply(RvF2,2,sum)
		RvCtr <- t( t(RvF2) / RvSF2)
	}
	plot.out <- prettyPlot(RvFS,constraints=constraints,col=participant.colors,main=ZeTitle,x_axis=axis1,y_axis=axis2,contributionCircles=TRUE,contributions=RvCtr,display_names=!nude)
	return(list(constraints=constraints,participant.colors=participant.colors))
}
