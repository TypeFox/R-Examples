`Amatdual` <-
function (steps, pointsin, removelist, nbrs, weights, alpha) 
{

n <- length(pointsin) + length(removelist)
Adual <- matrix(0, n - steps + 1, n - steps + 1)
Hdual <- matrix(0, n - steps, n - steps + 1)
Gdual <- matrix(0, 1, n - steps + 1)
newpoints <- (c(pointsin, rev(removelist)))[1:(n - steps + 1)]

o <-match(nbrs, newpoints[1:(length(newpoints) - 1)])
    
lastcol <- matrix(0, length(newpoints) - 1, 1)
lastcol[o] <- alpha
Hdual <- cbind(diag(length(newpoints) - 1), lastcol)

Hdual[o,o]<-as.column(-alpha)%*%as.row(weights)

diag(Hdual)[o]<-as.row(diag(Hdual)[o]+1)

Gdual[o] <- -weights
Gdual[length(newpoints)] <- 1
Adual <- rbind(Hdual, as.row(Gdual))

return(Adual)

}

