gpca <- function(X,Q,R,K=1,deflation = FALSE){

## error checking
	if(class(X) != "matrix" && class(X) != "dgCMatrix"){
		stop("X must be a matrix of class 'matrix'  or 'dgCMatrix' ")
	}
	
	if(class(Q) != "matrix" && class(Q) != "dgCMatrix"){
		stop("Q must be a matrix of class 'matrix'  or 'dgCMatrix' ")
	}
	
	
	if(class(R) != "matrix" && class(R) != "dgCMatrix"){
		stop("R must be a matrix of class 'matrix'  or 'dgCMatrix' ")
	}
	
	
	n = dim(X)[1]
	p = dim(X)[2]
	
	qn = dim(Q)[1]
	pn = dim(Q)[2]
	
	rn = dim(R)[1]
	rp = dim(R)[2]
	
	if(qn != pn){
		stop("Q must be a square matrix")
	}
	
	if(qn != n){
		stop("Q must be an n x n matrix")
	}
	
	if(rn != rp){
		stop("R must be a square matrix")
	}
	
	if(rn != p){
		stop("R must be a p x p matrix")
	}
		
	K = as.integer(K)
	if(K <1){
		stop("k must be an integer greater than 0")
	}
	
	return_list = 0
	
	return_list = tryCatch({
	
	if(deflation){
		return_list = gmdDeflation(X,Q,R,K,n,p)
	}else{
		
		if(n < p){
			temp = gmdLA(t(X),R,Q,K,p,n)
			return_list = list(temp[[2]],temp[[1]],temp[[3]],temp[[4]],temp[[5]])
		}else{
			return_list = gmdLA(X,Q,R,K,n,p)
		}
		
	}
	

			 
	},
			 warning = function(war){
				print(paste("MY_WARNING:  ",war))
				stop("Q and R must be Positive Semi Definite and have maximum eigen value less then or equal to 1")
			 },
			 error= function(err){
				print(paste("MY_ERROR:  ",err))
				stop("Q and R must be Positive Semi Definite and have maximum eigen value less then or equal to 1")
			 },
			 finally = {
			 
			 })
	
	names(return_list) = c("U","V","D","cumm.prop.var","prop.var")
	return(return_list)
	
}




