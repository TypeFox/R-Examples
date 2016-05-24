mywnn <-
function(train, test, weight){
	#implement general WNN with given weight for a new test; test must be a vector.
	
	train = as.matrix(train)
	n = dim(train)[1]
	d = dim(train)[2]-1
	X = as.matrix(train[, 1:d])
	Y = train[, d+1]
	Ysort = rep(0,n)
	
	dist = function(x){sqrt(t(x-test) %*% (x-test))}
	Dis = apply(X,1,dist)
	Ysort = Y[order(Dis)]
	
	label = ifelse(sum(weight[which(Ysort == 1)]) > 1/2, 1, 2)
	return(label)

}
