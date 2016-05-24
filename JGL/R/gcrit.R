gcrit <-
function(theta, S,n, lam1, lam2, penalize.diagonal)  # theta = list of pXp matrices, length k
{
	p = dim(S[[1]])[1]
	K = length(S)
	lam1 = penalty.as.matrix(lam1,p,penalize.diagonal=penalize.diagonal)
	lam2 = penalty.as.matrix(lam2,p,penalize.diagonal=penalize.diagonal)	
	crit = 0
	sumsquared = matrix(0,p,p)
	for(k in 1:length(theta))
	{
		# add log det that was entered as an argument, or else calculate it
		crit = crit+n[k]*log(det(theta[[k]]))-n[k]*sum(S[[k]]*theta[[k]])-sum(lam1*abs(theta[[k]])) 
		# add to term that will become group L2 norm:
		sumsquared = sumsquared + theta[[k]]^2
	}
	crit = crit - sum(lam2*sqrt(sumsquared))
	return(crit)
}

