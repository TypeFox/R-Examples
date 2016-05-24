"spec.emb" <- 
function(A,p,norm=TRUE){
		N <- nrow(A)
		if(all(A !=t(A))){
			A <- A + t(A)
		}
		DS <- diag(1/sqrt(apply(A,2,sum)))
		L <- diag(apply(A,1,sum)) -A
		NL <- DS%*%L%*%DS
		if(norm==TRUE){
			EIG <- eigen(NL)
		}else{
			EIG <- eigen(L)
		}
		return(EIG$vectors[,(N-p):(N-1)])
	}