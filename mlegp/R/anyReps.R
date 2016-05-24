`anyReps` <-
function(X) {
        X = as.matrix(X)
	for (i in 1:(dim(X)[1]-1)) {
		for (j in (i+1):dim(X)[1]) {
			if (identical(X[i,], X[j,])) return (TRUE)
		}
	}
	return (FALSE)
}

