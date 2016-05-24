#(size c.three way mixed cx(bina). model 5 a)
# Section 3.4.4.4 test factor A

# Three way mixed classification. (A>B)xC, Model V
# Factor A and B are fixed, factor C is random. 
# a, b and n are given, determining c. Testing hypothesis about factor A
size_c.three_way_mixed_cxbina.model_5_a <-  function(alpha, beta, delta, a, b, n, cases)
{
	c <- 5    
	c.new <- 1000
	while (abs(c -c.new)>1e-6)
	{
		c <- c.new
		dfn <- a-1
		dfd <- (a-1)*(c-1)
		lambda <- ncp(dfn,dfd,alpha,beta)
		if (cases == "maximin")
		{
			c.new <- 2*lambda/(b*n*delta*delta)
		}
		else if (cases == "minimin")
		{
			c.new <- 4*lambda/(a*b*n*delta*delta)
		}
	}  
	return(ceiling(c.new))
}


# example
# size.3_4_4_4.test_factor_A(0.05, 0.1, 0.5, 6, 5, 2, "maximin")
# size.3_4_4_4.test_factor_A(0.05, 0.1, 0.5, 6, 5, 2, "minimin")


