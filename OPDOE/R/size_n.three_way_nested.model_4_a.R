#(size n.three way nested. model 4 a)
# Section 3.4.2.3 test factor A

# Three way nested classification. Model IV
# Factor A and C fixed, B random. Determining n
# Testing hypothesis about factor A
size_n.three_way_nested.model_4_a <-  function(alpha, beta, delta, a, b, c, cases)
{
  if(is.na(b)){
    b <- 2
  n <- size_n.three_way_nested.model_4_a(alpha, beta, delta, a, b, c, cases)
  # print(b*n)
  min <- b*n
  b.min <- b
  n.min <- n
  while(!is.na(n)){
    b <- b+1 
    n <- size_n.three_way_nested.model_4_a(alpha, beta, delta, a, b, c, cases)
    if(is.na(n))break
    # print(n*b)
    if(n*b<min){
      min <- n*b
      n.min <- n
      b.min <- b
    }
  }
  cat(paste("\nminimum value of b*n is ",b.min*n.min,"\n\n"))
  list(b=b.min,n=n.min)
  } else {
	n <- 2
	dfn <- (a-1)
	dfd <- a*(b-1)
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
		n <- 5    
		n.new <- 1000
		while (abs(n -n.new)>1e-6)
		{
			n <- n.new
			dfn <- (a-1)
			dfd <- a*(b-1)
			lambda <- ncp(dfn,dfd,alpha,beta)
			if (cases == "maximin")
			{
				n.new <- 2*lambda/(b*c*delta*delta)
			}
			else if (cases == "minimin")
			{
				n.new <- 4*lambda/(a*b*c*delta*delta)
			}
		}  
		return(ceiling(n.new))
	}
      }
}


# example
# size.3_4_2_3.test_factor_A(0.05, 0.1, 0.5, 6, 5, 4, "maximin")
# size.3_4_2_3.test_factor_A(0.05, 0.1, 0.5, 6, 5, 4, "minimin")


optim_size_n.three_way_nested.model_4_a <- function(alpha, beta, delta, a, b, c, cases){
  n <- size_n.three_way_nested.model_4_a(alpha, beta, delta, a, b, c, cases)
  # print(b*n)
  min <- b*n
  b.min <- b
  n.min <- n
  while(!is.na(n)){
    b <- b+1 
    n <- size_n.three_way_nested.model_4_a(alpha, beta, delta, a, b, c, cases)
    if(is.na(n))break
    # print(n*b)
    if(n*b<min){
      min <- n*b
      n.min <- n
      b.min <- b
    }
  }
  cat(paste("\nminimum value of b*n is ",b.min*n.min,"\n\n"))
  list(b=b.min,n=n.min)
}
