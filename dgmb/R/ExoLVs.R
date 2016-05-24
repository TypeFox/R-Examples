ExoLVs <-
function(N,n,bv,bex,x.ex,weiex)
	{
	y.ex <- array(matrix(0,n,bex),dim=c(n,bex,N),dimnames=list(1:n,paste("y.ex",1:bex,sep=""),1:N))

	bv.acum <- cumsum(bv)
	for (i in 1:N)
		{
		aux <- 1
		for (j in 1:bex)
			{
			y.ex[,j,i] <- x.ex[,aux:bv.acum[j],i]%*%weiex[j,]
			aux <- bv.acum[j]+1
			}
		}
	list(y.ex=y.ex)
	}

