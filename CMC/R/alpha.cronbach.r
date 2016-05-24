alpha.cronbach = function(x){
	
	n = nrow(x) #n. of subject
	k = ncol(x) #n. of items

	var.cov.mat 	= cov(x) * (n-1)/n
	num 				= sum(diag(var.cov.mat))
	den 				= 2*sum(var.cov.mat[lower.tri(var.cov.mat)])+num
	alpha 			= k / (k-1) * (1 - num/den)
	return(alpha)
}