
vol.box <- function(box, x){
		
	x.ind.curr <- rep(TRUE, nrow(x))
  box.curr <- box

  for (j in 1:ncol(x)){
    x.ind.curr <- x.ind.curr & (x[,j]>= box.curr[1,j]) & (x[,j] <= box.curr[2,j])
	}

	return(sum(x.ind.curr))

}
		
#Original version:	
#`vol.box` <- function(box){
#  
#	#Note: 
#	return(prod(abs(box[2,] - box[1,])))
#
#}

##### Suggest this alternate: See github issue #15
#`vol.box` <- function(box){
#
#	#Make sure the box is well formed:
#	if(!all(box[1,] <= box[2,])){ 
#		cat("Problem detected in box formation: 
#		     one ore more box lower limits are higher than the higher limits or there are NAs.")
#		cat("Here is are the box limits, where the columns are the dimensions, the 
#		first row holds the lower limits and the second row holds the upper limits.")
#		print(box)
#		cat("The program will now give an error until that issue is dealt with.\n")
#		stop()
#	} 
#
#	return(prod(box[2,] - box[1,]))
#
#}

