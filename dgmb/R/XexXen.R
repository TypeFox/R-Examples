XexXen <-
function(N,n,ind.ex,ind.en,x.ex,x.en)
	{
		x <- array(matrix(0,n,ind.ex+ind.en),dim=c(n,ind.ex+ind.en,N))
		dimnames(x) <- list(1:n,c(paste("ind.ex",1:ind.ex,sep=""),paste("ind.en",1:ind.en,sep="")),1:N)
		x <- abind(x.ex,x.en,along=2)
		list(x=x)
	}

