data.blca <-
function(X){
	is.binary<- function(X) all(X[X>0]==1)	
	X<-as.matrix(X)
	M<- ncol(X)
	test1<- apply(X, 2, is.binary)
	
	if(sum(!test1)==1){
		Xnew<- NULL
		Xnew$counts.n<- X[, !test1]
		Xnew$data<- X[, test1]
		class(Xnew)<- "data.blca"
		return(Xnew)
	}
	
	if(sum(!test1)==0){
	  if(M>10){
	  	Xnew<- NULL
		Xnew$counts.n<- rep(1, nrow(X)) 
		Xnew$data<- X
		class(Xnew)<- "data.blca"
		return(Xnew)
	    }else{
	  N<- 2^M
	  counts.n<-countpattern(X)
			
	  if(!is.null(colnames(X))) colname1<- colnames(X)
	  else colname1<- NULL
			
	  X<-matrix(0, nrow=N, ncol=M)
	  colnames(X)<- colname1
			
	  for(i in 1:M) X[,M-i+1]<-rep(rep(c(0,1), c(2^(i-1), 2^(i-1))), 2^(M-i))
			
	  if(any(counts.n==0)) {
		X<-X[counts.n>0,]
		counts.n<-counts.n[counts.n>0] 
		}
			
		Xnew<- NULL
		Xnew$counts.n<- counts.n
		Xnew$data<- X
		class(Xnew)<- "data.blca"
		return(Xnew)
		}
		}else stop("Data is of incorrect type. It must be either an entirely  binary matrix, or a numeric matrix, with all columns binary, except for a 
			    single numeric column of counts.")
}
