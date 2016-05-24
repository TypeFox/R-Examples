
###############################################
# converts a nested list into a list
nestedList2List <- function( nestedList){
	NB <- length(nestedList)
	# count number of elements
	M <- 0
	for (bb in 1:NB){
	   M <- M + length(nestedList[[bb]])	
			}
	# create new list object
	list1 <- as.list(1:M)
	vv <- 1
	for (bb in 1:NB){
	  NW <- length(nestedList[[bb]])
	   for (ww in 1:NW){
			list1[[vv]] <- nestedList[[bb]][[ww]]
			vv <- vv + 1
			}
		}
	return(list1)	
	}
#################################################	