minmaxHelper <-
function(mat1,mat2=NULL,axis1=1,axis2=2,findBounds=TRUE){
	
	if(!is.matrix(mat1) && !is.data.frame(mat1)){
		warning(paste("mat1 was not a matrix or a data frame. mat1 class is: ",class(mat1),sep=""))
		return(0)
	}
	if(is.null(mat2)){
		mat2=mat1
	}

	if(findBounds){ #this gets the 2 mins & 2 maxs from x & y
		original_min_max_list <- findRealMinMax(mat1,mat2,axis1,axis2)
		minx<-original_min_max_list$minx * 1.15
		maxx<-original_min_max_list$maxx * 1.15	
		miny<-original_min_max_list$miny * 1.15	
		maxy<-original_min_max_list$maxy * 1.15
		
		if(minx == 0){
			minx = -abs(maxx * 0.1)
		}
		if(miny == 0){
			miny = -abs(maxy * 0.1)
		}	
		if(maxx == 0){
			maxx = abs(minx * 0.1)
		}		
		if(maxy == 0){
			maxy = abs(miny * 0.1)
		}
	}else{ #this finds the biggest (absolute) value for axis1 and axis2 in mat1 and mat2 and makes that the max/min
		max.val <- max(max(abs(mat1[,c(axis1,axis2)])),max(abs(mat2[,c(axis1,axis2)])))
		if(max.val == 0){
			#something must be wrong, but whatever.
			max.val <- 0.1
		}
		minx <- (-max.val * 1.15) -> miny
		maxx <- (max.val * 1.15) -> maxy		
	}
	minMaxList <- list(minx=minx,miny=miny,maxx=maxx,maxy=maxy)
	return(minMaxList)
}
