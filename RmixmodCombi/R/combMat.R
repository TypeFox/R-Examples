combMat <-
function(K,l1,l2)
{
	l=c(min(l1,l2), max(l1,l2))
	if(any(length(l1) == 0, length(l2) == 0)){
		l1 = numeric(0)
	  l2 = l[2]}
	 else {
	 	l1 = l[1]
	 	l2 = l[2]} 
	
	M <- rbind(cbind(diag(l2-1),matrix(rep(0,(K-l2+1)*(l2-1)), nrow=l2-1, ncol=K-l2+1)), cbind(matrix(rep(0,l2*(K-l2)), nrow=K-l2, ncol=l2),diag(K-l2)))
	M[l1,l2] <- 1
	return(M)
}
