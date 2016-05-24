#(size n.three way cross.model 1 a)
# Section 3.4.1.1 test factor A

# Three way cross classification - Model I
# Factor A, B, C are fixed. Determining sub-class number n
# Levels of factors A, B and C are given.
# Hypothesis: Factor A has no effects. 
size_n.three_way_cross.model_1_a <-  function(alpha, beta, delta, a, b, c, cases)
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
# size.3_4_1_1(0.05, 0.1, 0.5, 6, 5, 4, "maximin")
# size.3_4_1_1(0.05, 0.1, 0.5, 6, 5, 4, "minimin")


