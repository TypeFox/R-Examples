predict.smac = function(object,new.x=NULL,lambda=NULL,...)
{

if (is.null(new.x)) {new.x=object$x}

if (is.null(lambda)) {lambda=object$lambda}

if (class(new.x)!="matrix") {new.x = as.matrix(new.x)}
if (ncol(new.x)!=ncol(object$x)) {stop("The new covariates matrix/data.frame has a wrong dimension.")}
if (!is.numeric(lambda)) {stop("All lambda should be numeric.")}
if (any(lambda < 0)) {stop("All lambda should be non-negative.")}

pred.y=numeric(0)
pred.prob=numeric(0)
beta=numeric(0)
beta0=numeric(0)

	for (i in 1:length(lambda))	
	{
	temp=lambda[i]
	index=which(object$lambda==temp)
#############################################
	if (length(index)==1)
		{
		temp.beta=object$beta[[index]]
		temp.beta0=object$beta0[[index]]
		}
	if (length(index)==0)
		{
		if (temp>(object$lambda[1])) 
			{
			temp.beta=object$beta[[1]]
			temp.beta0=object$beta0[[1]]
			if (object$way==2) 
				{
				cat(paste("Lambda",temp,"is bigger than the largest lambda in the solution path.\nUsing the solution corresponding to the largest lambda in solution path.\n"))
				}
			}
		if (temp<(object$lambda[length(object$lambda)])) 
			{
			temp.beta=object$beta[[length(object$lambda)]]
			temp.beta0=object$beta0[[length(object$lambda)]]
			cat(paste("Lambda",temp,"is less than the smallest lambda in the solution path.\nUsing the solution corresponding to the smallest lambda in solution path.\n"))
			}
		if (temp>(object$lambda[length(object$lambda)]) & temp<(object$lambda[1]))
			{
			index2=max(which(object$lambda>temp))
			temp.beta=object$beta[[index2]]+(temp-object$lambda[index2])/(object$lambda[(index2+1)]-object$lambda[index2]) * (object$beta[[(index2+1)]]-object$beta[[index2]])
			temp.beta0=object$beta0[[index2]]+(temp-object$lambda[index2])/(object$lambda[(index2+1)]-object$lambda[index2]) * (object$beta0[[(index2+1)]]-object$beta0[[index2]])
			}
		}
#############################################

	temp.pred.y = pred.angle(new.x, temp.beta, temp.beta0, object$k)

	temp.pred.y2 = character(nrow(new.x))
	for (i in 1:object$k)
		{
		temp.pred.y2[temp.pred.y==i]=object$y.name[i]
		}

	if (is.numeric(object$y)) temp.pred.y2=as.numeric(temp.pred.y2)

	temp.pred.prob = prob.pred.angle(new.x, temp.beta, temp.beta0, object$k, object$loss)
	temp.pred.prob = data.frame(temp.pred.prob)
	colnames(temp.pred.prob) = object$y.name
	
	pred.y = c(pred.y,list(temp.pred.y2))
	pred.prob = c(pred.prob,list(temp.pred.prob))
	beta = c(beta,list(temp.beta))
	beta0 = c(beta0,list(temp.beta0))
	} #for (i in 1:length(lambda))	

a=list(new.x=new.x, lambda=lambda,fitted.beta0=beta0,fitted.beta=beta,pred.y=pred.y,pred.prob=pred.prob)

	this.call = match.call()
	a$call <- this.call

return(a)
}
