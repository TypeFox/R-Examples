#(size n.one way.model 1)
# section 3.2 one-way layout test factor A
#One way ANOVA: A fixed. determine n
size_n.one_way.model_1 <-  function(alpha, beta, delta, a, cases)
{
	n <- 2
        if(cases == "maximin" && a%%2==1) a=a-1
	dfn <- a-1
	dfd <- a*(n-1)
	if (cases == "maximin")
	{
		lambda <- 0.5*n*delta*delta
	}
	else if (cases == "minimin")
	{
		lambda <- 0.25*a*n*delta*delta
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
		n <- 5    
		n.new <- 1000
		while (abs(n -n.new)>1e-6)
		{
			n <- n.new
			dfn <- a-1
			dfd <- a*(n-1)
			lambda <- ncp(dfn,dfd,alpha,beta)
			if (cases == "maximin")
			{
				n.new <- 2*lambda/(delta*delta)
			}
			else if (cases == "minimin")
			{
				n.new <- 4*lambda/(a*delta*delta)
			}
		} 
		return(ceiling(n.new))
	}
}


# example
# size.one_way(0.05,0.1, 2, 4, "maximin")
# size.one_way(0.05,0.1, 2, 4, "minimin")


