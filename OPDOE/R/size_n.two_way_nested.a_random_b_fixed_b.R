#(size n.two way nested.a random b fixed b)
# Section 3.3.2.3 test factor B

# Two way nested classification A>B
# A random, B fixed. Testing hypothesis about factor B
# value of a is specified, determine n
size_n.two_way_nested.a_random_b_fixed_b <-  function(alpha, beta, delta, a, b, cases)
{
	n <- 2
	dfn <- a*(b-1)
	dfd <- a*b*(n-1)
	if (cases == "maximin")
	{
		lambda <- 0.5*n*delta*delta
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
			dfn <- a*(b-1)
			dfd <- a*b*(n-1)
			lambda <- ncp(dfn,dfd,alpha,beta)
			if (cases == "maximin")
			{
				n.new <- 2*lambda/(delta*delta)
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
# size.3_3_2_3(0.05, 0.1, 1, 2, 10, "maximin")
# size.3_3_2_3(0.05, 0.1, 1, 2, 10, "minimin")
# size.3_3_2_3(0.05, 0.1, 1, 3, 10, "maximin")
# size.3_3_2_3(0.05, 0.1, 1, 3, 10, "minimin")
# size.3_3_2_3(0.05, 0.1, 1, 10, 10, "maximin")
# size.3_3_2_3(0.05, 0.1, 1, 10, 10, "minimin")


