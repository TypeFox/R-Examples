#! Rank based Mahalanobis distance, from Paul Rosenbaum's Design of Observational Studies, p. 251
.smahal = function(z, X) { 
	X = as.matrix(X) 
	n = dim(X)[1] 
	rownames(X) = 1:n 
	k = dim(X)[2] 
	m = sum(z) 
	for (j in 1:k) X[,j] = rank(X[,j]) 
	cv = cov(X) 
	vuntied = var(1:n)
#! ***PENDING: correct this	
	diag(cv)[diag(cv) == 0] = .01
	rat = sqrt(vuntied/diag(cv)) 
	cv = diag(rat)%*%cv%*%diag(rat) 
	out = matrix(NA,m,n-m) 
	Xc = X[z == 0, ] 
	Xt = X[z == 1, ] 
	rownames(out) = rownames(X)[z==1] 
	colnames(out) = rownames(X)[z==0] 
	#library(MASS) 
	icov = ginv(cv) 
	for (i in 1:m) {
		out[i,] = mahalanobis(Xc,Xt[i,],icov,inverted=TRUE)
	} 
	out
}
