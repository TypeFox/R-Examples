"pexp" <-

function(q, rate = 1, beta = 1/rate,  lower.tail = TRUE, log.p = FALSE)



{



    if (any(beta <= 0) || any(rate <= 0)) 



            stop("beta (or rate) must be strictly positive")



	 stats:::pexp(q,rate=1/beta,lower.tail=lower.tail,log.p=log.p)

	

}

