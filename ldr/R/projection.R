projection <-
function(alpha, Sigma) 
{
	return(alpha%*%solve(t(alpha)%*%Sigma%*% alpha)%*%t(alpha)%*% Sigma)
}
