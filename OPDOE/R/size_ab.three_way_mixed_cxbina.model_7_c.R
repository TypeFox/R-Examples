#(size ab.three way mixed cx(bina). model 7 c)
size_ab.three_way_mixed_cxbina.model_7_c <-  function(alpha, beta, delta, c, n, cases)
{
    #assumption is: a>=2, b>=2, c>=2, n is fixed to 2
    #Note smaller b will result in larger a because both a and b are factors of lambda
    #Lambda is fixed for given precision
	
    ########### Step One ##########
    #Determine b.max letting a=2 
    a <- 2
    b.max <- size_b.three_way_mixed_cxbina.model_7_c(alpha, beta, delta, a, c, n, cases)
    print (b.max)
    resultN <- c(a,b.max,c,n)
    N.min <- 2*b.max*c*n
	
    ########### Step two ##########
    #for b from 2 to b.max determine a and calculate N.new=abcn, 
    #if N.new <N.old, memerize the value of a, b, c, and n
    for (b in 2:b.max)
    {
        a <- size_a.three_way_mixed_cxbina.model_7_c(alpha, beta, delta, b, c, n, cases)
        N.new <- a*b*c*n
        print (c(N.new, N.min, a, b, c, n))
        if (N.new < N.min)
        {
            resultN <- c(a,b,c,n)
            N.min <- N.new
        }
        if (a<2) 
        {
          break;
        }
    }
    return(ceiling(resultN))
}

# example
# size.3_4_4_7.find_minimalProduct(0.05,0.1,0.50, 5,2,  "maximin")
# size.3_4_4_7.find_minimalProduct(0.05,0.1,0.50, 5,2,  "minimin")



