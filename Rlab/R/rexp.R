"rexp" <-

function(n, rate = 1, beta = 1/rate)



{



    if (any(beta <= 0) || any(rate <= 0)) 



            stop("beta (or rate) must be strictly positive")



	 stats:::rexp(n,rate=1/beta)

	

}

