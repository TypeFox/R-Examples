#(size b.two way cross.mixed model a fixed a)
#  Section 3.3.1.2 test factor A

# Two way cross classification AxB
# A fixed, B random. determine b. 
# n can be fixed depending on whether there is an interaction
# If there is an interaction let n=2, otherwise let n=1
# Test equality of factor A
size_b.two_way_cross.mixed_model_a_fixed_a <-  function(alpha, beta, delta, a, n, cases)
{
	b <- 5    
	b.new <- 1000
	while (abs(b -b.new)>1e-6)
	{
		b <- b.new
		dfn <- a-1
		dfd <- (a-1)*(b-1)
		lambda <- ncp(dfn,dfd,alpha,beta)
		if (cases == "maximin")
		{
			b.new <- 2*lambda/(n*delta*delta)
		}
		else if (cases == "minimin")
		{
			b.new <- 4*lambda/(a*n*delta*delta)
		}
	}
	return(ceiling(b.new))
}


# example
# size.3_3_1_2(0.05,0.1, 1, 6, 4, "maximin")
# size.3_3_1_2(0.05,0.1, 1, 6, 4, "minimin")


