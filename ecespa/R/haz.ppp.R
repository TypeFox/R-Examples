`haz.ppp` <-
function (W) 
{
    
    if (dim(W)[2] == 2) 
        pepe = ppp(x = W[, 1], y = W[, 2], xrange = range(W[, 
            1]), yrange = range(W[, 2]))
    if (dim(W)[2] == 3) 
        pepe = ppp(x = W[, 1], y = W[, 2], xrange = range(W[, 
            1]), yrange = range(W[, 2]), marks = W[, 3])
    return(pepe)
}

