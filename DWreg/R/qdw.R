qdw<-function(p,q=exp(-1),beta=1)
{
if(q>1 | q <0)
	stop('q must be between 0 and 1', call. = FALSE)
if(beta <=0)
	stop('beta must be positive', call. = FALSE)
if(p<=1-q)
	res<-0
else
	res<-ceiling((log(1-p)/log(q))^(1/beta)-1)
return(res)
}