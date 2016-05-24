mybnn <-
function(train,test,ratio){
	#implement Bagged NN with given re-sampling ratio for new tests. test can be matrix or a vector
	# Reference: Hall and Samworth (2005)
	
	n = dim(train)[1]
	weight = rep(0,n)
	for(i in 1:n){
		weight[i] = ratio*(1-ratio)^(i-1)/(1-(1-ratio)^n)
	}
	
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
