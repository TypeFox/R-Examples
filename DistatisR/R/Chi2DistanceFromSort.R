Chi2DistanceFromSort <-
function(X){
	# iIve a Cube of Chi2Distance for a Sorting Task
	# X is a matrix of sort 
	# Each column is a participant
	# 2 objects in a column with the same number
	# were in the same group
	# send back a Object*Object*Participant
	# distance array
	
	# First Chi2Dist as a private function
	Chi2Dist    <-function(X){#compute the chi2 distance 
	# between the rows of a matrix
    # send back the distance and the vector of mass	
    
    # 1. transform X into Row profiles
    xip = rowSums(X)
    R <- X / xip
    # 2. Masses, weights
    xpp = sum(xip)             # grand total
    m <- xip / xpp             # masses 
    c <- colSums(X) / xpp      # row barycenter
    w = 1/c                    # columns weights
    # Preprocess R
    Rc = t(t(R) - c)           # deviations to barycenter
    Rtilde = t(t(Rc)*sqrt(w))  # weighted R
    S = Rtilde%*%t(Rtilde)     # covariance
    s =diag(S) # diag of
    D = (s - S) + t(s-S)       # Chi2 distance matrix
    return(list(Distance=D, masses=m))
} # end of private Chi2Dist

	nI = nrow(X);nJ = ncol(X)
	# initialize
	LeCube2Distance = array(0, c(nI,nI,nJ))
	# Horrible  Loop!
	for(j in 1:nJ){ 
		# get CHi2Distance of jth assessor expressed as disjunctive coding
		dist =  Chi2Dist(model.matrix(~ as.factor(as.matrix(X[,j])) - 1))
		LeCube2Distance[,,j] <- dist$Distance 
		 } # done ugly loop
	dimnames(LeCube2Distance) <-list(rownames(X),rownames(X),colnames(X))
	return(LeCube2Distance)  
}
