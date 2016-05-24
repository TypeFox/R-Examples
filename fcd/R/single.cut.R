single.cut <-
function(A, clusters, K = 2){
	
	ratio.count = 0
	normalised.count = 0
	temp.diag = apply(A, 2, sum)
	for(i in 1:K){
		temp.label = which(clusters == i)
		if(length(temp.label) != 0){
			t1 = sum(A[temp.label, -temp.label])
			ratio.count = ratio.count + t1/length(temp.label)
			normalised.count = normalised.count + t1/sum(temp.diag[temp.label])
		}
	}
	
	return(list(ratio.count = ratio.count, normalised.count = normalised.count))
}
