GraphDistatisCompromise <-
function(FS,axis1=1,axis2=2,constraints=NULL,item.colors=NULL,ZeTitle= 'Distatis-Compromise',nude=FALSE,Ctr=NULL){

    if(is.null(item.colors)){
		item.design <- diag(dim(FS)[1])
		item.colors <- as.matrix(createColorVectorsByDesign(item.design)$oc)
	}
    
	if(is.null(constraints) || sum(names(constraints) %in% c('minx','maxx','miny','maxy')) != 4){
		print('Making constraints')
    	#First, compute the constraints
	    real.minimum <- min(FS)
    	real.maximum <- max(FS)
	    real.value <- max(c(abs(real.minimum),abs(real.maximum)))
   		#set up a constraints list
		constraints <- list(minx=-real.value, maxx=real.value,miny=-real.value,maxy=real.value)
	}

	if(is.null(Ctr)){
		F2 <- FS**2
		SF2 <- apply(F2,2,sum)
		Ctr <- t( t(F2) / SF2)
	}
  # Now Compute the contributions for the plane of the dimensions to be plotted
  #Ctr4Plot =  apply(Ctr[,c(axis1,axis2)],1,sum)
  
  plot.out <- prettyPlot(FS,constraints=constraints,col=item.colors,main=ZeTitle,x_axis=axis1,y_axis=axis2,contributionCircles=TRUE,contributions=Ctr,display_names=!nude)

  return(list(constraints=constraints,item.colors=item.colors))
}
