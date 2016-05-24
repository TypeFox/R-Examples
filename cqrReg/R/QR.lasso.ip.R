QR.lasso.ip=function(X,y,tau,lambda)
{

if(missing(lambda)){
	lambda=1
}

a=rq.fit.lasso(X,y,tau,lambda)
return(a)
}
