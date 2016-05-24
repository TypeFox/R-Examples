### Enumerate wrapper function
## Input: Takes in margin vector MM and outcome vector NN
## Output: Matrix of possible internal cells for vector where N=1

enumerate = function(MM, NN){
	if(sum(MM) != sum(NN)){stop("MM and NN have different totals")}
	if(length(NN)!=2){stop("NN should be of length 2")}

      MMlen = length(MM)

	result <- .Call("Enumerate", MM, NN, MMlen, PACKAGE="osDesign")
	return(result)
}


enumerate.count = function(MM, NN){
	if(sum(MM) != sum(NN)){stop("MM and NN have different totals")}
	if(length(NN)!=2){stop("NN should be of length 2")}

      MMlen = length(MM)
	
	result <- .Call("EnumerateCount", MM, NN, MMlen, PACKAGE="osDesign")
	return(as.vector(result))
}

enumerate.window = function(MM, NN, startval, windowsize){
	if(sum(MM) != sum(NN)){stop("MM and NN have different totals")}
	if(length(NN)!=2){stop("NN should be of length 2")}

      MMlen = length(MM)

	result <- .Call("EnumerateWindow", MM, NN, MMlen, startval, windowsize, PACKAGE="osDesign")
	return(result)
}
