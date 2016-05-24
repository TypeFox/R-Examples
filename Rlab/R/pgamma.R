"pgamma" <-

function (q, shape = 1, rate = 1, scale = 1/rate, alpha = shape, beta = scale, lower.tail = TRUE, log.p = FALSE) 



{



    if (any(alpha <= 0)) 



            stop("alpha (or shape) must be strictly positive")

            

    if (any(beta <= 0) || any(scale <= 0) || any(rate <= 0)) 



            stop("beta (or scale or rate) must be strictly positive")



     stats:::pgamma(q,shape=alpha,scale=beta,lower.tail=lower.tail,log.p=log.p)



}

