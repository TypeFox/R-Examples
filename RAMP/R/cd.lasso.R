cd.lasso <-
function(X, y, a0, beta, epsilon, max.iter, lambda, family, bInd, nonPen){
	#dyn.load('cd_lasso1.so')
	###nonPen p-dimensional indicator whether it impose a penalty or not. nonPen=1: no penalty
	X = cbind(X,1)
	p = ncol(X)
	n = nrow(X)
	beta = c(beta,a0)
	nonPenAll = c(nonPen,rep(0,p-1-length(nonPen)),1)
	#nonPenAll[1:q] = nonPen
	if(family=='binomial') family=2
	if(family=='poisson') family=1
	if(family=='gaussian') family=0
	para.in = c(epsilon, max.iter, lambda, family)
	#cat('beta[1:10]', beta[1:min(length(beta),10)],'\n')
	if(family == 0){
	out = .Fortran('cd_lasso_lin',
	       X = as.double(X),
	       y = as.double(y),
	       p = as.integer(p),
	       n = as.integer(n),
	       beta = as.double(beta),
	       nonPen = as.integer(nonPenAll),
	       paraIn = as.double(para.in)     
	       )      	
	       }
	else if(family == 2) {
	out = .Fortran('cd_lasso_bin',
	       X = as.double(X),
	       y = as.double(y),
	       p = as.integer(p),
	       n = as.integer(n),
	       beta = as.double(beta),
	       nonPen = as.integer(nonPenAll),
	       paraIn = as.double(para.in)     
	       )            	
	       	
	       }       
return(list(a0=out$beta[p], beta=out$beta[-p]))		
}
