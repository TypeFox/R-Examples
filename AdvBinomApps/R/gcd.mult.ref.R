gcd.mult.ref <- function(A, prec=2){
	
	r<-length(A)
	A.prec<-round(A,prec)	
	A.prec[A.prec-A>0]<-A.prec[A.prec-A>0]-10^(-prec)
	A<-round(A.prec*10^prec)
	A<-sort(A,decreasing=TRUE)

	gcd<-function(a,b){
		if(b==0){
			return<-a
		}else{
			return<-gcd(b,a%%b)
		}
		return
	}

	A.gcd<-A[1]
	if(r>1){
		for(i in c(2:r)){
			A.gcd<-gcd(A.gcd,A[i])
		}
	}
	A.gcd/10^prec
}
