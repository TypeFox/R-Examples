#(size a.three way mixed cx(bina). model 3 c)
# Section 3.4.4.2 test factor C

# Three way mixed classification. (A>B)xC, Model III
# Factor A is random, factor B and C are fixed. n = 1.
# b and c are given, determining a. Testing hypothesis about factor C
size_a.three_way_mixed_cxbina.model_3_c <-  function(alpha, beta, delta, b, c, n, cases)
{
	a <- 5    
	a.new <- 1000
	while (abs(a - a.new)>1e-6)
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
	}  
	return(ceiling(a.new))
}


# example
# size.3_4_4_2.test_factor_C(0.05, 0.1, 0.5, 5, 4, 1, "maximin")
# size.3_4_4_2.test_factor_C(0.05, 0.1, 0.5, 5, 4, 1, "minimin")



