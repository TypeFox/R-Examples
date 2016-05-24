GraphDistatisBoot <-
function(FS,FBoot, axis1=1, axis2=2, item.colors=NULL, ZeTitle= 'Distatis-Bootstrap', constraints=NULL, nude = FALSE, Ctr=NULL, lwd = 3.5, ellipses=TRUE, fill = TRUE, fill.alpha = .27, percentage=0.95){         
	
	if (is.null(item.colors)){
		item.design <- diag(dim(FS)[1])
		item.colors <- as.matrix(createColorVectorsByDesign(item.design)$oc)
	}

	if(is.null(constraints)){
		real.minimum <- min(c(FS,FBoot))
    	real.maximum <- max(c(FS,FBoot))
	    real.value <- max(c(abs(real.minimum),abs(real.maximum)))
   		#set up a constraints list
		constraints <- list(minx=-real.value, maxx=real.value,miny=-real.value,maxy=real.value)
	}

	Compromise.Out <- GraphDistatisCompromise(FS,axis1=axis1,axis2=axis2,constraints=constraints,item.colors=item.colors,ZeTitle=ZeTitle,nude=nude,Ctr=Ctr)
    for(i in 1:dim(FS)[1]){    
    	if(ellipses){	
	        dataEllipse(x=FBoot[i,axis1,],y=FBoot[i,axis2,],add=TRUE,levels=percentage,plot.points=FALSE,center.pch = FALSE, col = item.colors[i,], lwd = lwd, fill=fill,fill.alpha=fill.alpha)	
		}else{
			##this intentionally overlays a black background hull.
			peeledHull(t(FBoot[i,,]),x_axis=axis1,y_axis=axis2,percentage=percentage,lwd=lwd)
			if(lwd<3){
				color.lwd <- 1
			}else{
				color.lwd <- lwd-2
			}
			peeledHull(t(FBoot[i,,]),x_axis=axis1,y_axis=axis2,percentage=percentage,col=item.colors[i,],lwd=color.lwd)	
		}
	}
    return(Compromise.Out)   
#return(LeRetour)
}
