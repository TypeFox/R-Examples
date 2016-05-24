rlca <-
function(n, itemprob=0.5, classprob=1, fit=NULL){ 

if(is.null(fit))
{
	itemprob<- as.matrix(itemprob)
	G<- nrow(itemprob)
	M<- ncol(itemprob)
	} else{
		itemprob<- fit$itemprob
		classprob<- fit$classprob
		G<- nrow(itemprob)
		M<- ncol(itemprob)
	}

x <- matrix(runif(n*M), nrow=n)
classvec<- as.vector(rmultinom(1, n, prob=classprob))

ind<- c(0, cumsum(classvec))

for(g in 1:G) x[(ind[g]+1):ind[g+1],] <- t(t(x[(ind[g]+1):ind[g+1],]) < itemprob[g,])*1

return(x)

}
