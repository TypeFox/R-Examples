EnLVs <-
function(N,n,ben,bet,elv,y.ex.tot)
	{
		bet <- as.matrix(bet)	
		y.en <- array(matrix(0,n,ben),dim=c(n,ben,N),dimnames=list(1:n,paste("y.en",1:ben,sep=""),1:N))							
		for (i in 1:N)  y.en[,ben,i] <- y.ex.tot[,,i]%*%bet+elv[,,i]
		list(y.en=y.en)
	}

