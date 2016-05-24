#(size b.three way mixed ab in c. model 3 a)
# Section 3.4.3.2 test factor A

# Three way mixed classification. (AxB)>C, Model III
# Factor A, C are fixed. B is random. Determining b, n=1
# a and c are given. Testing hypothesis about factor A
size_b.three_way_mixed_ab_in_c.model_3_a <-  function(alpha, beta, delta, a, c, n, cases)
{
	b <- 5    
	b.new <- 1000
	while (abs(b -b.new)>1e-6)
	{
		b <- b.new
		dfn <- (a-1)
		dfd <- (a-1)*(b-1)
		lambda <- ncp(dfn,dfd,alpha,beta)
		if (cases == "maximin")
		{
			b.new <- 2*lambda/(c*n*delta*delta)
		}
		else if (cases == "minimin")
		{
			b.new <- 4*lambda/(a*c*n*delta*delta)
		}
	}  
	return(ceiling(b.new))
}


# example
# size.3_4_3_2.test_factor_A(0.05, 0.1, 0.5, 6, 5, 1, "maximin")
# size.3_4_3_2.test_factor_A(0.05, 0.1, 0.5, 6, 5, 1, "minimin")



