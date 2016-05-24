`varPerReps` <-
function(X, Y) {
	X = as.matrix(X)
	Y = as.matrix(Y)
	if (!anyReps(X)) stop("no reps detected in estimateNugget")

  	complete = matrix(FALSE, ncol=1, nrow = dim(X)[1])

	i = 1
	v = matrix(0,ncol=1, nrow = dim(X)[1])
  	for (i in 1:dim(X)[1]) {
		if (complete[i]) next

        	index = matrix(FALSE, dim(X)[1])
		for (j in 1:dim(X)[1]) {
			index[j] = all(X[i,] == X[j,])
		}
		if (sum(index) > 1) {
			v[index] = var(Y[index])
		}
		else {
			v[index] = 0
		}
      		complete = complete | index
	}
	return (v)
}

