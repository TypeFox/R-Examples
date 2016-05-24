mysnn <-
function(train, test, lambda){
	#implement SNN with lambda by predicting output for a new test
	
	n = dim(train)[1]
	d = dim(train)[2]-1
	X = as.matrix(train[, 1:d])
	Y = train[, d+1]
	Ysort = rep(0,n)
	weightstar = rep(0,n)
	
	Kstar = floor(((d*(d+4)/(2*(d+2)))*lambda*n^(4/d))^(d/(d+4)))
	if(Kstar == 0){Kstar = 1}else if(Kstar > n){ Kstar = n}
	
	for(i in 1:Kstar){
		weightstar[i] = (1+d/2-d/(2*Kstar^(2/d))*(i^(1+2/d)-(i-1)^(1+2/d)))/Kstar
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
	
	label = rep(0,nrow(test.mat))
	
	for(j in 1:nrow(test.mat)){
	
		dist=function(x){sqrt(t(x-test.mat[j,])%*%(x-test.mat[j,]))}
		Dis=apply(X,1,dist)
		Ysort = Y[order(Dis)]
		label[j] = ifelse(sum(weightstar[which(Ysort[1:Kstar]==1)])>1/2,1,2)
	}
	
	return(label)

}
