near3 <-
function (x, data) 
{#find the near row to x in data using manhattan distance
    nd <- length(data[, 1])
    distall <- rep(0, nd)
    for (i in 1:nd) {
        distall[i] <- distancia1(x, data[i, ])
    }
    ind1 <- order(distall)[1]
    near1 <- data[ind1, ]
    near1
}
