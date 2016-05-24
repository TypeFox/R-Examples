
# Tests whether two sets of shifts are exactly identical
#	a = vector of shift nodes
#   b = vector of shift nodes

areShiftSetsEqual <- function(a, b){
	if (length(a) != length(b)){
		return(FALSE);
	}else if (length(a) == 0 & length(b) == 0){
		return(TRUE);
	}else{
		if (length(intersect(a,b)) != length(a)){
			return(FALSE);
		}else{
			return(TRUE);
		}
	}
}

 


