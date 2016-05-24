#(size b.three way mixed cx(bina). model 4 c)
# Section 3.4.4.3 test factor C

# Three way mixed classification. (A>B)xC, Model IV
# Factor A and C are fixed, factor B is random. n = 1.
# a and c are given, determining b. Testing hypothesis about factor C
size_b.three_way_mixed_cxbina.model_4_c <-  function(alpha, beta, delta, a, c, n, cases)
{
	b <- 5    
	b.new <- 1000
	while (abs(b -b.new)>1e-6)
	{
		b <- b.new
		dfn <- c-1
		dfd <- a*(b-1)*(c-1)
		lambda <- ncp(dfn,dfd,alpha,beta)
		if (cases == "maximin")
		{
			b.new <- 2*lambda/(a*n*delta*delta)
		}
		else if (cases == "minimin")
		{
			b.new <- 4*lambda/(a*c*n*delta*delta)
		}
	}
	return(ceiling(b.new))
}


# example
# size.3_4_4_3.test_factor_C(0.05, 0.1, 0.5, 6, 4, 1, "maximin") 
# size.3_4_4_3.test_factor_C(0.05, 0.1, 0.5, 6, 4, 1, "minimin")


