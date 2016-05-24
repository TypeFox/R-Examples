ErrEnLV <-
function(N,n,ben,bet,y.ex.cor)
	{
	std.error <- rep(0,N)
	vis <- 0
	for (i in 1:N)
		{
		aux <- 1-t(bet)%*%y.ex.cor[,,i]%*%bet
		ifelse(aux < 0, vis <- 1, std.error[i] <- sqrt(1-t(bet)%*%y.ex.cor[,,i]%*%bet))
		}
	if (vis != 0) print("warnings NAs, the sample size is too small, please set a higher sample size")
	if (vis != 0) stop
	
	elv <- array(matrix(0,n,ben),dim=c(n,ben,N),dimnames=list(1:n,paste("d",1:ben,sep=""),1:N))						
	for (i in 1:N)
		{
		elv[,,i] <- apply(as.matrix(elv[,,i]), 2, function(x) rnorm(n,0,std.error[i]))
		}

	list(std.error=std.error,elv=elv,vis=vis)
	}

