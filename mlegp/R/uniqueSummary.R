uniqueSummary <-function (X, Y) 
{
    X = as.matrix(X)
    Y = as.matrix(Y)
    if (!anyReps(X)) { 
        #stop("no reps detected in meanPerReps")
        cat("no reps detected in meanPerReps, returning originals...\n")
	reps = rep(1, nrow(X))
	uniqueX = X
	uniqueMeans = Y
	uniqueVar = rep(0,nrow(X))
	l = list(reps = reps, uniqueX = uniqueX, uniqueMeans = uniqueMeans, uniqueVar = uniqueVar)
  	return(l)
	}

   uniqueX = unique(X)
   uniqueMeans = NULL
   uniqueVar = NULL
   reps = NULL
   for (i in 1:nrow(uniqueX)){
        index = matrix(FALSE, dim(X)[1])
        for (j in 1:dim(X)[1]) {
            index[j] = all(uniqueX[i, ] == X[j, ])
        }
	
		uniqueMeans = rbind(uniqueMeans, apply(matrix(Y[index,],ncol=ncol(Y)), 2, mean))
		uniqueVar =   rbind(uniqueVar,   apply(matrix(Y[index,],ncol=ncol(Y)), 2, var))

	reps = rbind(reps, sum(index))
   }

  l = list(reps = reps, uniqueX = uniqueX, uniqueMeans = uniqueMeans, uniqueVar = uniqueVar)
  return(l)
}
