#(size c.three way mixed cx(bina). model 7 b)
# Section 3.4.4.6 test factor B in A

# Three way mixed classification. (A>B)xC, Model VII
# Factor A and C are random, factor B is fixed. 
# a, b and n are given, determining c. Testing hypothesis about factor B in A
size_c.three_way_mixed_cxbina.model_7_b <-  function(alpha, beta, delta, a, b, n, cases)
{
	c <- 2
	dfn <- a*(b-1)
	dfd <- a*(b-1)*(c-1)
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
	        c <- 5    
	        c.new <- 1000
	        while (abs(c -c.new)>1e-6)
	        {
		c <- c.new
		dfn <- a*(b-1)
		dfd <- a*(b-1)*(c-1)
		lambda <- ncp(dfn,dfd,alpha,beta)
		if (cases == "maximin")
		{
			c.new <- 2*lambda/(n*delta*delta)
		}
		else if (cases == "minimin")
		{
			c.new <- 4*lambda/(a*b*n*delta*delta)
		}
	}  
	return(ceiling(c.new))
        }
}


# example
# size.3_4_4_6.test_factor_BinA(0.05, 0.1, 0.5, 6, 5, 2, "maximin")
# size.3_4_4_6.test_factor_BinA(0.05, 0.1, 0.5, 6, 5, 2, "minimin")



#######################################################################







#######################################################################

