DistanceFromSort <-
function(X){
	# Private functions first
	# Sort2Dist: Create a sorting distance matrix from the sorting vector
    Sort2Dist <- function(LeSort){# Start Sort2Dist
	 nObj = length(LeSort)
     truc = matrix(as.matrix(LeSort),nrow=nObj,ncol=nObj)
     # in lieu of repmat!
     DistMat = 1 - (truc==t(truc)) # to get 0/1 instead of T/F
	 return(DistMat)
     }# End Sort2Dist
	# X is a matrix of sort 
	# Each column is a participant
	# 2 objects in a column with the same number
	# were in the same group
	# send back a Object*Object*Participant
	# distance array
	nI = nrow(X);nJ = ncol(X)
	# initialize
	LeCube2Distance = array(0, c(nI,nI,nJ))
	for(j in 1:nJ){LeCube2Distance[,,j]<-Sort2Dist(X[,j]) }
	dimnames(LeCube2Distance) <-list(rownames(X),rownames(X),colnames(X))
	return(LeCube2Distance)  
}
