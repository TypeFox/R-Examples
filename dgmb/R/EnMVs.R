EnMVs <-
function(N,n,ind.en,weien,y.en)
	{
	x.en <- array(matrix(0,n,ind.en),dim=c(n,ind.en,N),dimnames=list(1:n,paste("x.en",1:ind.en,sep=""),1:N))

	for (i in 1:N)
		{
		x.en[,2:ind.en,i] <- apply(as.matrix(x.en[,2:ind.en,i]),2,function(x) rnorm(n,0,1))
		}

	for (i in 1:N)
		{
		x.en[,1,i] <- (1/weien[1,1])*(y.en[,,i]-as.matrix(x.en[,2:ind.en,i])%*%weien[1,2:ind.en])
		}
	list(x.en=x.en)
	}

