
###########################################
# convert a List into a nestedList
List2nestedList <- function( List , N_between , N_within = NULL ,
		loop_within = TRUE ){

	M <- length(List)	
	if ( is.null(N_within)){ N_within <- floor(M / N_between) }
	
	nl <- as.list(1:N_between)
	M2 <- N_between * N_within
	for (bb in 1:N_between){
		nl[[bb]] <- as.list(1:N_within)
			}
	
	
	#****************************
	# loop within
	bb <- 1
	ww <- 1
	if ( loop_within){
	for (mm in 1:M2){		
		nl[[bb]][[ww]] <- List[[ mm ]]		
		if ( ww < N_within){ 
				ww <- ww + 1
					} else {
				bb <- bb + 1 
				ww <- 1
					}		
			}
		}
	#****************************
	# loop between
	if ( ! loop_within){
	for (mm in 1:M2){		
		nl[[bb]][[ww]] <- List[[ mm ]]		
		if ( bb < N_between){ 
				bb <- bb + 1
					} else {
				ww <- ww + 1 
				bb <- 1
					}		
			}
		}		
	return(nl)
		}
#########################################################		