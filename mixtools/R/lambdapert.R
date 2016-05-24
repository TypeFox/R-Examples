lambda.pert <- function (lambda, pert) 
{
	temp = logit(lambda) + pert
	temp2 = inv.logit(temp)
	new.lambda = temp2/sum(temp2)
	new.lambda
}
