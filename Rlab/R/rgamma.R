"rgamma" <-

function (n, shape = 1, rate = 1, scale = 1/rate, alpha = shape, beta = scale) 



{



    if (any(alpha <= 0)) 



            stop("alpha (or shape) must be strictly positive")

            

    if (any(beta <= 0) || any(scale <= 0) || any(rate <= 0)) 



            stop("beta (or scale or rate) must be strictly positive")



     stats:::rgamma(n,shape=alpha,scale=beta)



}

