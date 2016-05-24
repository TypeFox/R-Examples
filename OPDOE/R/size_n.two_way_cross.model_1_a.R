#(size n.two way cross.model 1 a)
# Section 3.3.1.1  test factor A

# Two way cross classification AxB
# A and B fixed. determine n
# Testing equality of factor A
size_n.two_way_cross.model_1_a <-  function(alpha, beta, delta, a, b, cases)
{
	n <- 2
	dfn <- a-1
	dfd <- a*b*(n-1)
	if (cases == "maximin")
	{
		lambda <- 0.5*b*n*delta*delta
	}
	else if (cases == "minimin")
	{
		lambda <- 0.25*a*b*n*delta*delta
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
		n <- 5    
		n.new <- 1000
		while (abs(n -n.new)>1e-6)
		{
			n <- n.new
			dfn <- a-1
			dfd <- a*b*(n-1)
			lambda <- ncp(dfn,dfd,alpha,beta)
			if (cases == "maximin")
			{
				n.new <- 2*lambda/(b*delta*delta)
			}
			else if (cases == "minimin")
			{
				n.new <- 4*lambda/(a*b*delta*delta)
			}
		}  
		return(ceiling(n.new))
	}
}


# example
# size.3_3_1_1(0.05,0.1, 1, 6, 4, "maximin")
# size.3_3_1_1(0.05,0.1, 1, 6, 4, "minimin")


