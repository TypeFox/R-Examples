#(size n.three way mixed cx(bina). model 1 b)
# Section 3.4.4.1 test factor B

# Three way mixed classification. (A>B)xC, Model I
# Factor A, B and C are fixed. Determining n.
# a b and c are given. Testing hypothesis about factor B in A
size_n.three_way_mixed_cxbina.model_1_b <-  function(alpha, beta, delta, a, b, c, cases)
{
	n <- 2
	dfn <- a*(b-1)
	dfd <- a*b*c*(n-1)
	if (cases == "maximin")
	{
		lambda <- 0.5*c*n*delta*delta
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
		n <- 5    
		n.new <- 1000
		while (abs(n -n.new)>1e-6)
		{
			n <- n.new
			dfn <- a*(b-1)
			dfd <- a*b*c*(n-1)
			lambda <- ncp(dfn,dfd,alpha,beta)
			if (cases == "maximin")
			{
				n.new <- 2*lambda/(c*delta*delta)
			}
			else if (cases == "minimin")
			{
				n.new <- 4*lambda/(a*b*c*delta*delta)
			}
		}
		return(ceiling(n.new))
	}
}


# example
# size.3_4_4_1.test_factor_B(0.05, 0.1, 0.5, 6, 5, 4, "maximin")
# size.3_4_4_1.test_factor_B(0.05, 0.1, 0.5, 6, 5, 4, "minimin")

