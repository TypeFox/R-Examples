Inv.proj.matrix.sdp <-
function(matrix){
	mat.decomp			      	=	eigen(matrix,symmetric=TRUE)
	valpp.mat			        	=	mat.decomp$values
	valpp.mat[which(valpp.mat<0.00001)]	=	0
	valpp.mat_inv						=	sapply(valpp.mat,function(X){ifelse(X==0,0,1/X)})
	if(length(valpp.mat_inv)==1){
		valpp.mat_inv=as.matrix(valpp.mat_inv)
	}
	mat.decomp$vectors %*% (diag(valpp.mat_inv)) %*% t(mat.decomp$vectors)
}
