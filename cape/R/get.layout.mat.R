get.layout.mat <-
function(num.panes, type = c("landscape", "upright")){
	
	#first find the square root of the number of panes
	edge.dim <- sqrt(num.panes)
	round.dim <- round(edge.dim)
	
	#if the number is square, just return a square matrix
	if(round.dim == edge.dim){
		layout.mat <- matrix(1:num.panes, ncol = edge.dim, byrow = TRUE)
		return(layout.mat)
		}
	
	#if the the edge.dim is not a whole number
	# we need to do some more figuring
	#if the nearest whole number is less than
	#the square root of the number of panes
	#add one to get close
	if(round.dim < edge.dim){
		dims <- c(round.dim, (round.dim+1))
		}else{
			dims <- c(round.dim, (round.dim - 1)) #otherwise, subtract one to get close
			}

		total.panes <- prod(dims)
		test.dims <- matrix(rep(dims, 2), nrow = 2, byrow = TRUE)			
		
		#as long as we haven't found a spot for all plots
		while(all(total.panes < num.panes)){
			#add rows and columns until we get
			#equal to or greater than the total
			#number of panes.
			for(i in 1:length(test.dims[,1])){
				max.dim.locale <- which(test.dims[1,] == max(test.dims[1,]))
				test.dims[i,max.dim.locale[1]] <- test.dims[i,max.dim.locale[1]] + 1
				}		
			total.panes <- apply(test.dims, 1, prod)
			}
			
			
			#figure out if either of these two options equals the original num.panes
			found.it <- which(total.panes == num.panes)
			if(length(found.it) > 0){
				found.it.locale <- which(total.panes == num.panes)
				dims <- test.dims[found.it.locale[1],]
				}else{
					#otherwise, figure out which of these is the minimum of 
					#the values that goes over the total number of panes
					greater <- which(total.panes > num.panes)
					dim.locale <- which(total.panes == min(total.panes[greater]))
					dims <- test.dims[dim.locale[1],]
					}
			
			
			#Make the layout matrices
			#fill in 0's for the panes we don't
			#want to fill with figures
			
			#the default layout is landscape
			if(length(grep("l", type)) > 0){
				layout.mat <- matrix(c(1:num.panes, rep(0, (total.panes[1]-num.panes))), ncol = max(dims), byrow = TRUE)
				}else{
				layout.mat <- matrix(c(1:num.panes, rep(0, (total.panes[1]-num.panes))), ncol = min(dims), byrow = TRUE)
				}
		
	return(layout.mat)
	
	
	
	
}
