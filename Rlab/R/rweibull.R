"rweibull" <-

function (n, shape = 1, scale = 1, alpha = shape, beta = scale) 



{



    if (any(alpha <= 0)) 



            stop("alpha (or shape) must be strictly positive")

            

    if (any(beta <= 0)) 



            stop("beta (or scale) must be strictly positive")



     stats:::rweibull(n,shape=alpha,scale=beta)



}

