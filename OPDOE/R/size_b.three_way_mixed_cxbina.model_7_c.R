#(size b.three way mixed cx(bina). model 7 c)
size_b.three_way_mixed_cxbina.model_7_c <-  function(alpha, beta, delta, a, c, n, cases)
{
	b <- 2
	dfn <- (c-1)
	dfd <- (a-1)*(c-1)
	if (cases == "maximin")
	{
		lambda <- 0.5*a*b*n*delta*delta
	}
	else if (cases == "minimin")
	{
		lambda <- 0.25*a*b*c*n*delta*delta
	}
	beta.calculated <- Beta(alpha, dfn, dfd, lambda)
	if (is.nan(beta.calculated) || beta.calculated < beta )
	{
   warning(paste("Given parameter will result in too high power.",
                 "To continue either increase the precision or ",
                 "decrease the level of factors."))
               return(NA)
	}
	else
	{
		b <- 2   
		b.new <- 1000
		count.loop <- 1
		while (abs(b -b.new)>1e-2 && count.loop<20)
		{
			b <- b.new
			dfn <- (c-1)
			dfd <- (a-1)*(c-1)
			lambda <- ncp(dfn,dfd,alpha,beta)
			if (cases == "maximin")
			{
				b.new <- 2*lambda/(a*n*delta*delta)
			}
			else if (cases == "minimin")
			{
				b.new <- 4*lambda/(a*c*n*delta*delta)
			}
			count.loop <- count.loop + 1
		}
		return(b.new)
	}
}


