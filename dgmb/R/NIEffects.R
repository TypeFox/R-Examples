NIEffects <-
function(N,n,y.ex,rie)
	{
	a.nle <- array(matrix(0,n,ncol(rie)),dim=c(n,ncol(rie),N),dimnames=list(1:n,paste("nle",1:ncol(rie),sep=""),1:N))					
	a.ie <- array(matrix(0,n,ncol(rie)),dim=c(n,ncol(rie),N),dimnames=list(1:n,paste("ie",1:ncol(rie),sep=""),1:N))

	for (i in 1:N)										
		{
		for (k in 1:nrow(rie))
			{
			for (j in 1:ncol(rie))
				{
				if (rie[k,j]==1 & k == j) a.nle[,k,i] <- y.ex[,k,i]*y.ex[,k,i]
				}
			}
		}
	for (i in 1:N)										
		{
		aux <- 1
		for (k in 1:nrow(rie))
			{
			for (j in 1:ncol(rie))
				{
				if (rie[k,j]==1 & k != j)
					{
					a.ie[,aux,i] <- y.ex[,k,i]*y.ex[,j,i]
					aux <- aux + 1
					}
				}
			}
		}
	list(a.nle=a.nle,a.ie=a.ie)
	}

