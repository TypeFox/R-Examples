#(size b=c.three way mixed cx(bina). model 6 a case2)
# Section 3.4.4.5 test factor A, Sigma(b(a))=0 , b=c, case 2

# Three way mixed classification. (A>B)xC, Model VI
# Factor A is fixed, factor B and C are random. 
# a and n are given, determining b and c. Testing hypothesis about factor A
# case2: assumption Sigma(b(a)) = 0 , b=c
size_bc.three_way_mixed_cxbina.model_6_a_case2 <-  function(alpha, beta, delta, a, n, cases)
{
	b <- 2
	c <- b
	dfn <-  a-1
	dfd <- (a-1)*(c-1)
	if (cases == "maximin")
	{
		lambda <- 0.5*b*c*n*delta*delta
	}
	else if (cases == "minimin")
	{
		lambda <- 0.25*a*b*c*n*delta*delta
	}
	beta.calculated <- Beta(alpha, dfn, dfd, lambda)
	if (is.nan(beta.calculated) || beta.calculated < beta )
	{
   warning(paste("Given parameter will result in too high power.",
                 "To continue either increase the precision or ",
                 "decrease the level of factors."))
               return(NA)
	}
	else
	{
	        b <- 5    
	        b.new <- 1000
                c <- b
	        while (abs(b -b.new)>1e-6)
	        {
		b <- b.new
                c <- b
		dfn <-  a-1
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
}


# example
# size.3_4_4_5.test_factor_A.case2(0.05, 0.1, 0.5, 6,  2, "maximin")
# size.3_4_4_5.test_factor_A.case2(0.05, 0.1, 0.5, 6,  2, "minimin")


