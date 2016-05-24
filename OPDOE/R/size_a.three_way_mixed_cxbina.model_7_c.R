#(size a.three way mixed cx(bina). model 7 c)
# Section 3.4.4.7 test factor C, implementation which find the minimal N=abcn

# Three way mixed classification. (A>B)xC, Model VIII
# Factor A and B are random, factor C is fixed.
# c is given, n=2, determining b for given a.
# Testing hypothesis about factor C
size_a.three_way_mixed_cxbina.model_7_c <-  function(alpha, beta, delta, b, c, n, cases)
{
	a <- 2
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
		a <- 2   
		a.new <- 1000
		count.loop <- 1
		while (abs(a -a.new)>1e-2 && count.loop<20)
		{
			a <- a.new
			dfn <- (c-1)
			dfd <- (a-1)*(c-1)
			lambda <- ncp(dfn,dfd,alpha,beta)
			if (cases == "maximin")
			{
				a.new <- 2*lambda/(b*n*delta*delta)
			}
			else if (cases == "minimin")
			{
				a.new <- 4*lambda/(b*c*n*delta*delta)
			}
			count.loop <- count.loop + 1
		} 
		return(a.new)
	}
}

