myknn <-
function(train, test, K){
	#implement KNN with given K for a new test; test can be matrix or a vector
	
	n = dim(train)[1]
	weight = rep(0,n)
	weight[1:K] = 1/K
	
    if(is.vector(test) == TRUE){
        
        if(dim(train)[2] - 1 == 1){
            # d = 1 case
            test.mat = as.matrix(test)
            
        }else{
            # d > 1 case
            test.mat = t(as.matrix(test))
        }
		
	}else{
	
		test.mat = test
	}
	
	if(dim(test.mat)[2] != (dim(train)[2]-1)) stop("training data and test data have different dimensions")	
	
	label = apply(test.mat,1,function(x) mywnn(train,x,weight))
	return(label)

}
