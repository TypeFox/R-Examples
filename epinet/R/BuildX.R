#######################
# Function  BuildX
#   Inputs: 
#			nodecov, an N x k matrix where each row represent and individual, column 1 is the id and the remaining (k-1) columns are covariate values for the node
#			unaryCol, an array of column indices
#			unaryFunc, an array of same length as unaryCol containing methods for comparing dyads.  Each entry can take values in: "match", "absdiff"
#			binaryCol, a list of 2 element arrays of column indices
#			binaryFunc, an array of same length as binaryCol containing methods for comparing dyads.  Each entry can take values in: "euclidean", "manhattan"
#			includeIntercept a logical value.  If TRUE, includes a column of all ones.  Defaults to TRUE
#   Outputs:
#			a dyadic covariate matrix with (N choose 2) rows, columns 1 and 2 are node ids, column 3 is all ones (if requested) and then one column for each given element of unaryCol and binaryCol.  
#			Builds colnames depending on type of unaryFunc and binaryFunc and colnames of nodecov
# example:
# 			mycov = data.frame(id = 1:5,xpos = rnorm(5),ypos = rnorm(5),house = c(1,1,2,2,2),gender = c(0,0,0,1,1))
# 			dyadCov = BuildX(mycov,unaryCol = c(4,5),unaryFunc = c("match","match"),binaryCol = list(c(2,3)),binaryFunc = c("euclidean"))
########################


BuildX <- function(nodecov,unaryCol = NULL,unaryFunc = NULL, binaryCol = NULL,binaryFunc = NULL,includeIntercept = TRUE){
	# check there at least one column
	if (!(includeIntercept | (length(unaryCol)>0) | (length(binaryCol)>0) )) {
		stop("need to have at least one column in X Matrix (either includeIntercept must be true or one of unaryCol or binaryCol must be non-empty)")
	}
  
	# check vector lengths match up
	if (length(unaryCol) != length(unaryFunc)){
		stop("unaryCol and unaryFunc must have same length")
	}
	if (length(binaryCol) != length(binaryFunc)){
		stop("binaryCol and binaryFunc must have same length")
	}
	
  # helper function
  comb2 <- function(x){
    ## given x with length(x)=n, returns array with n choose 2 rows and 2 columns giving all possible pairs of x
    ## faster version of combn(x,2)
    n  = length(x)
    k = n * (n-1) / 2
    y = array(rep(1:(n-1), (n-1):1), dim = c(k,2))
    for (i in 1:(n-1)) {
      z = k - ((n-i+1)*(n-i)/2)
      y[(z + 1):(z + n - i),2] =  (i+1):n
    }
    return(array(x[y], dim = c(n*(n-1)/2, 2)))
  }
  
	# n is number of nodes
	n = nrow(nodecov)
	covname = colnames(nodecov)
	# ncol is number of columns in final matrix
	ncol = 2 + includeIntercept + length(unaryCol) + length(binaryCol)
	
	# initialise final matrix
	dyadmat = matrix(0,nrow = n/2 *(n-1),ncol = ncol)
	colname = as.character(1:ncol)
	
	# make all dyads
	# dyadmat[,1:2] = t(combn(nodecov[,1],2))
	dyadmat[,1:2] = comb2(nodecov[,1])
	colname[1:2] = c("node.1","node.2")
	
	#keep track of how many columns filled
	completedCol = 2
	
	if (includeIntercept){
		dyadmat[,3] = 1
		colname[3] = "(Intercept)"
		completedCol = 3
	}
	
	if (length(unaryCol) > 0) {
		for (i in 1:length(unaryCol)) {
			currCol = nodecov[,unaryCol[i]]
			if (unaryFunc[i] == "match"){
				dyadmat[,completedCol+1] = currCol[dyadmat[,1]] == currCol[dyadmat[,2]]
				colname[completedCol+1] = paste(covname[unaryCol[i]],"match",sep = ".")
			} else if (unaryFunc[i] == "absdiff"){
				dyadmat[,completedCol+1] = abs(currCol[dyadmat[,1]] - currCol[dyadmat[,2]])
				colname[completedCol+1] = paste(covname[unaryCol[i]],"diff",sep = ".")
			} else {
				stop(paste("Entry",i,",",unaryFunc[i],", in unaryFunc not recognised"))
			}
			completedCol = completedCol + 1;
		}
	}

	if (length(binaryCol) > 0) {
		for (i in 1:length(binaryCol)) {
			currCol = nodecov[,binaryCol[[i]]]
			if (binaryFunc[i] == "euclidean"){
				dyadmat[,completedCol+1] = sqrt((currCol[dyadmat[,1],1] - currCol[dyadmat[,2],1])^2 + (currCol[dyadmat[,1],2] - currCol[dyadmat[,2],2])^2)
				colname[completedCol+1] = paste(covname[binaryCol[[i]][1]],covname[binaryCol[[i]][2]],"L2Dist",sep = ".")
			} else if (binaryFunc[i] == "manhattan"){
				dyadmat[,completedCol+1] = abs(currCol[dyadmat[,1],1] - currCol[dyadmat[,2],1]) + abs(currCol[dyadmat[,1],2] - currCol[dyadmat[,2],2])
				colname[completedCol+1] = paste(covname[binaryCol[[i]][1]],covname[binaryCol[[i]][2]],"L1Dist",sep = ".")
			} else {
				stop(paste("Entry",i,",",binaryFunc[i],", in binaryFunc not recognised"))
			}
			completedCol = completedCol + 1;
		}
	}
	
	colnames(dyadmat) <- colname

	return(dyadmat)
}