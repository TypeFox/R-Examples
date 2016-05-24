#(size c.three way cross.model 3 a)
# Section 3.4.1.2 test factor A

# Three way cross classification - Model III
# Factor A and B fixed, C random. n is fixed to 2. 
# a and b are given. Determining c
# Hypothesis: Factor A has no effects. 
size_c.three_way_cross.model_3_a <-  function(alpha, beta, delta, a, b, n, cases)
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
# size.3_4_1_2(0.05, 0.1, 0.5, 6, 5, 2, "maximin")
# size.3_4_1_2(0.05, 0.1, 0.5, 6, 5, 2, "minimin")


