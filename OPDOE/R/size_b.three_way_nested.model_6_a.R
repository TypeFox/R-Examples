#(size b.three way nested. model 6 a)
# Section 3.4.2.5 test factor A

# Three way nested classification. Model VI
# Factor A fixed, B and C random. Determining b,
# a and c are given. Testing hypothesis about factor A
size_b.three_way_nested.model_6_a <-  function(alpha, beta, delta, a, c, n, cases)
{
    b <- 2
    dfn <- a-1
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
         b <- 2   
         b.new <- 1000
         while (abs(b -b.new)>1e-5)
         {
                b <- 0.5*(b+ b.new)
                dfn <- (a-1)
                dfd <- a*(b-1)
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
}


# example
# size.3_4_2_5.test_factor_A(0.05, 0.1, 0.5, 6, 4, 2, "maximin")
# size.3_4_2_5.test_factor_A(0.05, 0.1, 0.5, 6, 4, 2, "minimin")


