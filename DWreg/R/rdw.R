rdw<-function(n,q=exp(-1),beta=1)
{
if(q>1 | q <0)
	stop('q must be between 0 and 1', call. = FALSE)
if(beta <=0)
	stop('beta must be positive', call. = FALSE)
y<-runif(n,0,1)
x<-as.numeric(lapply(y,qdw,q=q,beta=beta))
return(x)
}