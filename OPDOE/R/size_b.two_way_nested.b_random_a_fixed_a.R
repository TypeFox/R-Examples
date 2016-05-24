#(size b.two way nested.b random a fixed a)
# Section 3.3.2.2 test factor A

# Two way nested classification A>B
# A fixed, B random. Test factor A. Determining level b of Factor B
# In this model n=1
size_b.two_way_nested.b_random_a_fixed_a <-  function(alpha, beta, delta, a, cases)
{
	b <- 5    
	b.new <- 1000
	while (abs(b -b.new)>1e-6)
	{
		b <- b.new
		dfn <- a-1
		dfd <- a*(b-1)
		lambda <- ncp(dfn,dfd,alpha,beta)
		if (cases == "maximin")
		{
			b.new <- 2*lambda/(delta*delta)
		}
		else if (cases == "minimin")
		{
			b.new <- 4*lambda/(a*delta*delta)
		}
	}
	return(ceiling(b.new))
}


# example
# size.3_3_2_2(0.05, 0.1, 1, 6, "maximin")
# size.3_3_2_2(0.05, 0.1, 1, 6, "minimin")



