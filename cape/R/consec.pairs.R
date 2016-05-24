consec.pairs <-
function(elements){
	
	pair.mat <- matrix(NA, nrow = (length(elements)-1), ncol = 2)
	pair.mat[,1] <- elements[1:(length(elements)-1)]
	pair.mat[,2] <- elements[2:length(elements)]
	return(pair.mat)
	}
