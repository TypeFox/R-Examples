ProgressBAR<-function(n,N,...){
	n=n/N
	for(i in c(0, n)){
		A<-txtProgressBar(...)	
		setTxtProgressBar(A,i)
	}
}