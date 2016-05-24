percentile95 <-
function(X){		

n=ncol(X)
Q025=vector("numeric",length=n)
Q975=vector("numeric",length=n)
for (i in 1:n){
	Q025[i]=quantile(X[,i],0.025)
	Q975[i]=quantile(X[,i],0.975)
}

out=list()
out$lo=Q025
out$up=Q975
return(out)
}
