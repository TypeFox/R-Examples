ExoMVs <-
function(N,n,ind.ex)
	{
	x.ex <- array(matrix(0,n,ind.ex),dim=c(n,ind.ex,N),dimnames=list(1:n,paste("x.ex",1:ind.ex,sep=""),1:N))
	for (i in 1:N)
		{
		x.ex[,,i] <- apply(x.ex[,,i],2,function(x) rnorm(n,0,1))
		}
	list(x.ex=x.ex)
}

