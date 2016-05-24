fstep.GramSc <-
function(Y,T,kernel){

	n = nrow(Y)
	p = ncol(Y)
	K = ncol(T)
	m = colMeans(Y)
	d = min(p-1,(K-1))
	U = matrix(NA,p,d)

	# Compute S
	XX = as.matrix(Y - t(m*t(matrix(1,n,p))))
	TT = t(apply(T,1,"/",sqrt(colSums(T))))
	
	if (n>p & kernel==''){
		S = cov(Y)*(n-1)/n
		B = t(XX)%*%TT%*%t(TT)%*%XX / n
	
		# Eigendecomposition of S^-1 * B
		eig = eigen(solve(S)%*%B)
		U[,1] = matrix(Re(eig$vec[,1]),ncol=1)

if (d>1) for (k in 2:d){ 
			W = as.matrix(U[,-c(k:d)])  
			base = diag(p)
			base <- as.matrix(base[,1:(p-k+1)])
			W =  cbind(W,base)
	
	v <- W  

		# start of the gram-schmidt algorithm
		for (l in 2:p){
				proj <- c(crossprod(W[,l],v[,1:(l-1)])) / diag(crossprod(v[,1:(l-1)]))
				v[,l] <- matrix(W[,l],ncol=1) - matrix(v[,1:(l-1)],nrow=p) %*%  matrix(proj,ncol=1)
				v[,l] <- v[,l] /  sqrt(crossprod(v[,l]))
		}


	P   = v[,k:ncol(v)] 
	B.p = crossprod(P,crossprod(t(B),P))
	S.p = crossprod(P,crossprod(t(S),P))
	eig = eigen(ginv(S.p)%*%B.p)

	Umax = matrix(Re(eig$vec[,1]),ncol=1) 
        U[,k]= P%*%Umax	
	}
}
if (d==1) {U = U[,1]} 
as.matrix(U)
}

