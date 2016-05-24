`negloglik.laplace` <-
function (xpar, xx) 
{
    a <- xpar[2]
    w <- wfromt(xpar[1], a = a)
    loglik <- sum(log(1 + w * beta.laplace(xx, a)))
	if (abs(loglik)==Inf){
		loglik<-sign(loglik)*10^6
	}   
 return(-loglik)
}

