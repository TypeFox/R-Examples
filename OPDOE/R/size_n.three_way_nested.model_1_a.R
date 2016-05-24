#(size n.three way nested. model 1 a)
# Section 3.4.2.1 test factor A

# Three way nested classification. Model I
# Factor A, B and C are fixed. Determining n
# Testing hypothesis about factor A
size_n.three_way_nested.model_1_a <-  function(alpha, beta, delta, a, b, c, cases)
{
	n <- 2
	dfn <- a-1
	dfd <- a*b*c*(n-1)
	if (cases == "maximin")
	{
		lambda <- 0.5*b*c*n*delta*delta
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
			dfn <- a-1
			dfd <- a*b*c*(n-1)
			lambda <- ncp(dfn,dfd,alpha,beta)
			if (cases == "maximin")
			{
				n.new <- 2*lambda/(b*c*delta*delta)
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
# size.3_4_2_1.test_factorA(0.05, 0.1, 0.5, 6, 5, 4, "maximin")
# size.3_4_2_1.test_factorA(0.05, 0.1, 0.5, 6, 5, 4, "minimin")


