isolate <- function(A){
	n = dim(A)[1]
	temp = apply(abs(A), 2, sum)
	D = diag(temp, nrow = n)
	temp1 = which(apply(D, 2, sum) == 0)
	temp2 = which(apply(D, 2, sum) != 0)
	return(list(isolate = temp1, nonisolate = temp2))
}

