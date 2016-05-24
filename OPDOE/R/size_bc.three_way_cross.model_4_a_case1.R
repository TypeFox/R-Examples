#(size b=c.three way cross.model 4 a case1)
# Section 3.4.1.3 test factor A, sigma(AC)=0, case 1

#Three way cross classification AxBxC
#A is fixed. B and R are random
#Assumption is SigmaAC=0. b=c, n=2
#Testing effect of factor A
#Determining b=c
size_bc.three_way_cross.model_4_a_case1 <-  function(alpha, beta, delta, a, n, cases)
{
	b <- 5    
	b.new <- 1000
	while (abs(b -b.new)>1e-6)
	{
		b <- b.new
		c <- b
		dfn <- a-1
		dfd <- (a-1)*(b-1)
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
# size.3.4.1.3.case1(0.05, 0.1, 0.5, 6, 2, "maximin")
# size.3.4.1.3.case1(0.05, 0.1, 0.5, 6, 2, "minimin")


