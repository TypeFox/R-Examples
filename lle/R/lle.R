lle <-
function(X,m,k,reg=2,ss=FALSE,p=0.5,id=FALSE,nnk=TRUE,eps=1,iLLE=FALSE,v=0.99) {
	
	#vector of chosen data (subset selection)
	choise <- c() #remains empty if there's no subset selection
	
	#find neighbours
	if( iLLE ) cat("finding neighbours using iLLE\n") else cat("finding neighbours\n") 
	if( nnk ) nns <- find_nn_k(X,k,iLLE)	else nns <- find_nn_eps(X,eps)
		
	#calculate weights
	cat("calculating weights\n")
	res_wgts <- find_weights(nns,X,m,reg,ss,p,id,v)
	
	#if subset selection, the neighbours and weights have to be 
	#calculated again since the dataset changed
	if( ss ){
		
		#use reduced dataset
		X <- res_wgts$X
		choise <- res_wgts$choise
		
		cat("finding neighbours again\n")
		if( nnk ) nns <- find_nn_k(X,k,iLLE)	else nns <- find_nn_eps(X,eps)
		
		cat("calculating weights again\n")
		res_wgts <- find_weights(nns,X,m,reg,FALSE,p,id,v)
	}
	
	#extracting data (weights, intrnsic dim)
	wgts <- res_wgts$wgts
	id <- res_wgts$id
	
	#compute coordinates
	cat("computing coordinates\n")
	Y <- find_coords(wgts,nns,N=dim(X)[1],n=dim(X)[2],m)
	 
	return( list(Y=Y,X=X,choise=choise,id=id) )
}

