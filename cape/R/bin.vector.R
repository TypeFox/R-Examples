bin.vector <-
function(vectorX, bins = seq(0,1,0.5)){
	
	dist.mat <- apply(matrix(bins, ncol = 1), 1, function(x) x - vectorX)	
	binned.vector <- apply(dist.mat, 1, function(x) bins[which(abs(x) == min(abs(x)))[1]])
	return(binned.vector)	
	
}
