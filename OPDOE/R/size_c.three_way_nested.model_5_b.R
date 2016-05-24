#(size c.three way nested. model 5 b)
# Section 3.4.2.4 test factor B in A

# Three way nested classification. Model V
# Factor A and B fixed, C random. Determining c, n is fixed to 1
# Testing hypothesis about factor B in A
size_c.three_way_nested.model_5_b <-  function(alpha, beta, delta, a, b, n, cases)
{
	c <- 5    
	c.new <- 1000
	while (abs(c -c.new)>1e-6)
	{
		c <- c.new
		dfn <- a*(b-1)
		dfd <- a*b*(c-1)
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


# example
# size.3_4_2_4.test_factor_BinA(0.05, 0.1, 0.5, 6, 5, 1, "maximin")
# size.3_4_2_4.test_factor_BinA(0.05, 0.1, 0.5, 6, 5, 1, "minimin")



