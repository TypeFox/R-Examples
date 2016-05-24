"pweibull" <-

function (q, shape = 1, scale = 1, alpha = shape, beta = scale, lower.tail=TRUE, log.p = FALSE) 



{



    if (any(alpha <= 0)) 



            stop("alpha (or shape) must be strictly positive")

            

    if (any(beta <= 0)) 



            stop("beta (or scale) must be strictly positive")



     stats:::pweibull(q,shape=alpha,scale=beta,lower.tail=lower.tail,log.p=log.p)



}

