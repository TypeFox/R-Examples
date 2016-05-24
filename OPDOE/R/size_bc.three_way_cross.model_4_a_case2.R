#(size b=c.three way cross.model 4 a case2)
# Section 3.4.1.3 factor A, sigma(AB) = 0, case 2

#Three way cross classification AxBxC
#A is fixed. B and R are random
#Assumption is SigmaAB=0. b=c, n=2
#Testing effect of factor A
#Determining b=c
size_bc.three_way_cross.model_4_a_case2 <-  function(alpha, beta, delta, a, n, cases)
{
	b <- 5    
	b.new <- 1000
	c <- b
	while (abs(b -b.new)>1e-6)
	{
		b <- b.new
		c <- b
		dfn <- a-1
		dfd <- (a-1)*(c-1)
		lambda <- ncp(dfn,dfd,alpha,beta)
		if (cases == "maximin")
		{
			b.new <- sqrt(2*lambda/(n*delta*delta))
		}
		else if (cases == "minimin")
		{
			b.new <- sqrt(4*lambda/(a*n*delta*delta))
		}
	}
	return(ceiling(b.new))
}


# example
# size.3.4.1.3.case2(0.05, 0.1, 0.5, 6, 2, "maximin")
# size.3.4.1.3.case2(0.05, 0.1, 0.5, 6, 2, "minimin")


