myownn <-
function (train, test, K) {
	#implement S12's ownn with test K as that from KNN.
	# Reference: Samworth (2012), AOS.
	
	n = dim(train)[1]
	d = dim(train)[2]-1

    kstar <- floor((2 * (d + 4)/(d + 2))^(d/(d + 4)) * K)
	if(kstar > n){kstar = n}
    i <- 1:kstar
    alpha <- i^(1 + 2/d) - (i - 1)^(1 + 2/d)
    weight = c((1 + d/2 - d/(2 * kstar^(2/d)) * alpha)/kstar, rep(0, n - kstar))
	
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
