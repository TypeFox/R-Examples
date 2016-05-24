gmdLA <- function(X,Q,R,k,n,p){
	
##computation
	
	decomp = eigen(R)
	Rtilde = decomp$vectors %*% diag(sqrt(decomp$values)) %*% t(decomp$vectors)
	inv.values = decomp$values
	inv.values[which(decomp$values!=0)] = 1/sqrt(decomp$values[which(decomp$values!=0)])
	Rtilde.inv = decomp$vectors %*% diag(inv.values) %*% t(decomp$vectors)
	inmat =  t(X) %*% Q %*% X
	
	RtinRt = Rtilde %*% inmat %*% Rtilde
	xtilde.decomp = eigen(RtinRt)
	
	vgmd = Rtilde.inv %*% xtilde.decomp$vectors
	dgmd = sqrt(xtilde.decomp$values[1:k])
	
	ugmd = matrix(nrow = n,ncol = k)
	cumv = rep(0,k)
	
	
	propv = dgmd^2/sum(diag(as.matrix(inmat %*% R)))
	
	normalizing.number = 1
	
	XR = X %*% R
	RnR = R %*% inmat %*% R
	
	for( i in 1:k){
		normalizing.number = sqrt(vgmd[,i] %*% RnR %*% vgmd[,i])
		ugmd[,i] = as.matrix(XR %*% vgmd[,i])/as.double(normalizing.number)
		cumv[i] = sum(propv[1:i]) 
	}
		
	return.object = list("u"=ugmd,"v"=vgmd[,1:k],"d"=dgmd,"cumv"=cumv,"propv"=propv)
	return(return.object)

}




