myglmnet.max <-
function(X, link ="gaussian",delta=0.01){
	minlambda =0;
	#maxlambda = lambdaMax(t(X));
	
	tmp = X %*% t(X)
	maxlambda = max(tmp[upper.tri(tmp)])
	
	
	### binary search the interval
	while(1){
		mid = (minlambda+maxlambda)/2
		#tmp=glmpois(X,mid)
		tmp=glm.network(X, NULL, link=link, lambda=mid)		
		tmp[abs(tmp)<1e-06]=0
		tmp[abs(tmp)>1e-06]=1
		
		if(sum(tmp)>0){
			minlambda = mid+delta
		}else{
			maxlambda = mid-delta
		}			
		
		if(abs(maxlambda-minlambda)<delta){
			return(mid);
		}
	}
}
