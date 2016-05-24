
################################################
# multiple grouping helper functions
BIFIE_create_pseudogroup <- function( datalistM , 
	group , group_index , group_values ){

	GR <- length(group)
	res0 <- list( "datalistM" = datalistM , 
				"group_orig" = group , 
				"group"=group , "group_index"=group_index ,				
				"GR"= GR , "group_values" = group_values)
						
				
	#****************	
	#*** multiple groupings
	if (GR>1){
	    group_values <- as.list( 1:GR )
			for (gg in 1:GR){
				t1 <- fasttable( datalistM[ , group_index[gg] ] )				  
				group_values[[gg]] <- sort( as.numeric( paste( names(t1) ) ))
							}	
		res0$group_values_orig <- group_values
		
		datalistM2 <- datalistM[ , group_index]
		for (gg in 1:GR){
			datalistM2[,gg] <- match( datalistM2[,gg] , group_values[[gg]] )
						}
		maxval_exp <- 3				
		maxval <- 10^maxval_exp			
		res0$maxval <- maxval
		pseudogroup <- datalistM2[,1]
		for (gg in 2:GR){
			pseudogroup <- pseudogroup + maxval^(gg-1) * datalistM2[,gg]
						}
		t1 <- fasttable( pseudogroup )				  
		group_values <- sort( as.numeric( paste( names(t1) ) ))		
		res0$group_values <- group_values
				
		#**** group values recalculated in original values
		group_values_recode <- matrix( NA , nrow=length(group_values) , ncol=GR )
			
		for (gg in 1:GR){
			group_values_recode[,gg] <- group_values / maxval^(GR-gg) 			
						}		
		for (gg in 1:GR){
			group_values_recode[,gg] <- round( group_values_recode[,gg] , 0 )
					}			
		for (gg in 2:GR){
				group_values_recode[,gg] <- group_values_recode[,gg] %% ( maxval ) 			
					}
		group_values_recode <- group_values_recode[ , seq(GR,1 ,-1) ]		
		for (gg in 1:GR){
			h1 <- res0$group_values_orig[[gg]]
			group_values_recode[,gg] <- h1[ group_values_recode[,gg] ]			
						}
		res0$group_values_recode <- group_values_recode		
		res0$datalistM <- as.matrix( cbind( datalistM , pseudogroup ) )
		res0$group_index <- ncol(datalistM)+1
		res0$group <- "pseudogroup"
		
			}
	#****************


	return(res0)
	}
###################################################


