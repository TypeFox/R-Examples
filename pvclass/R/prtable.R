prtable <-
function(Y, pv, alpha = 0.05)
{
	pr <- (pv > alpha)
	L <- dim(pr)[2]
	if (L == 2)
	{
		T <- matrix(0,2,4)
		if (!is.null(dimnames(pv)[[2]]))
		{
			dimnames(T) <- list(b=dimnames(pv)[[2]],
				c("P(b,{})","P(b,{1})","P(b,{2})","P(b,{1,2})"))
		}
		else
		{
			dimnames(T) <- list(b=1:2,
				c("P(b,{})","P(b,{1})","P(b,{2})","P(b,{1,2})"))
		}
		for (b in 1:2)
		{
			ib <- (Y == b)
			T[b,1] <- sum(!pr[ib,1] & !pr[ib,2])
			T[b,2] <- sum( pr[ib,1] & !pr[ib,2])
			T[b,3] <- sum(!pr[ib,1] &  pr[ib,2])
			T[b,4] <- sum( pr[ib,1] &  pr[ib,2])
			T[b,]  <- T[b,] / sum(ib)
		}
	}
	if (L == 3)
	{
		T <- matrix(0,3,8)
		if (!is.null(dimnames(pv)[[2]]))
		{
			dimnames(T) <- list(b=dimnames(pv)[[2]],
				c("P(b,{})","P(b,{1})","P(b,{2})","P(b,{3})",
					"P(b,{1,2})","P(b,{1,3})","P(b,{2,3})",
					"P(b,{1,2,3})"))
		}
		else
		{
			dimnames(T) <- list(b=1:3,
				c("P(b,{})","P(b,{1})","P(b,{2})","P(b,{3})",
					"P(b,{1,2})","P(b,{1,3})","P(b,{2,3})",
					"P(b,{1,2,3})"))
		}
		for (b in 1:3)
		{
			ib <- (Y == b)
			T[b,1] <- sum(!pr[ib,1] & !pr[ib,2] & !pr[ib,3])
			T[b,2] <- sum( pr[ib,1] & !pr[ib,2] & !pr[ib,3])
			T[b,3] <- sum(!pr[ib,1] &  pr[ib,2] & !pr[ib,3])
			T[b,4] <- sum(!pr[ib,1] & !pr[ib,2] &  pr[ib,3])
			T[b,5] <- sum( pr[ib,1] &  pr[ib,2] & !pr[ib,3])
			T[b,6] <- sum( pr[ib,1] & !pr[ib,2] &  pr[ib,3])
			T[b,7] <- sum(!pr[ib,1] &  pr[ib,2] &  pr[ib,3])
			T[b,8] <- sum( pr[ib,1] &  pr[ib,2] &  pr[ib,3])
			T[b,]  <- T[b,] / sum(ib)
		}
	}
  if (L > 3)
  {
  		T <- matrix(0,L,L+3)
  		txt <- NULL
  		for (b in 1:L)
  		{
  			txt <- c(txt,paste("I(b,",as.character(b),")",sep=""))
  		}
  		txt <- c(txt,"P(b,{b})","P(b,{})",
  			paste("P(b,{1,...,",as.character(L),"})",sep=""))
		if (!is.null(dimnames(pv)[[2]]))
		{
			dimnames(T) <- list(b=dimnames(pv)[[2]],txt)
		}
		else
		{
			dimnames(T) <- list(b=1:L,txt)
		}
		for (b in 1:L)
		{
			ib <- (Y == b)
			for (theta in 1:L)
			{
				T[b,theta] <- sum(pr[ib,theta])
			}
			T[b,L+1] <- sum(pr[ib,b] &
				(apply(pr[ib,],1,sum) == 1))
			T[b,L+2] <- sum((apply(pr[ib,],1,sum) == 0))
			T[b,L+3] <- sum((apply(pr[ib,],1,sum) == L))
			T[b,]  <- T[b,] / sum(ib)
		}
  	}

  cat('\n')
  print(T)
  cat('\n')
  return(T)
}
