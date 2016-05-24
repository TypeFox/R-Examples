dsgl <-
function(A,L,lam1,lam2,penalize.diagonal)
{
	lam1 = lam1*1/L
	lam2 = lam2*1/L

  if(is.matrix(A[[1]])) {p=dim(A[[1]])[1]}
  if(is.vector(A[[1]])) {p=length(A[[1]])}
	K=length(A)
	softA = A
	for(k in 1:K) {softA[[k]] = soft(A[[k]],lam1,penalize.diagonal=penalize.diagonal) }   #if penalize.diagonal=FALSE was used in ggl(), then this will not penalize the diagonal.
	normsoftA = A[[1]]*0
	for(k in 1:K) {normsoftA = normsoftA + (softA[[k]])^2}

	normsoftA = sqrt(normsoftA)

	notshrunk = (normsoftA>lam2)*1
	# reset 0 elements of normsoftA to 1 so we don't get NAs later. 
	normsoftA = normsoftA + (1-notshrunk)

	out = A
	for(k in 1:K)
	{
		out[[k]] = softA[[k]]*(1-lam2/normsoftA)
		out[[k]] = out[[k]]*notshrunk
	}
	return(out)
}

