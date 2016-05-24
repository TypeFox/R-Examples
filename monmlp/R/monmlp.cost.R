monmlp.cost <-
function(weights, x, y, hidden1, hidden2, Th, To, Th.prime, To.prime,
         monotone=NULL)
{
    w <- monmlp.reshape(x=x, y=y, weights=weights, hidden1=hidden1,
                        hidden2=hidden2)
    W1 <- w$W1; W2 <- w$W2
    if (hidden2 > 0) W3 <- w$W3
    if (!is.null(monotone)){
        W1[monotone,] <- exp(W1[monotone,])
        W2[1:(nrow(W2)-1),] <- exp(W2[1:(nrow(W2)-1),])
        if(hidden2 > 0) W3[1:(nrow(W3)-1),] <- exp(W3[1:(nrow(W3)-1),])
    }
    x <- cbind(x, 1)
    h1 <- x %*% W1
    y1 <- Th(h1)
    aug.y1 <- cbind(y1, 1)
    h2 <- aug.y1 %*% W2
    if (hidden2==0){
        y2 <- To(h2)
        E <- y-y2
        deltaE2 <- 2*To.prime(h2)*E
        g.W2 <- -(t(aug.y1) %*% deltaE2)/length(E)
        if(!is.null(monotone))
            g.W2[1:(nrow(W2)-1),] <- g.W2[1:(nrow(W2)-1),]*W2[1:(nrow(W2)-1),]
        E1 <- deltaE2 %*% t(W2[1:(nrow(W2)-1),,drop=FALSE])
        deltaE1 <- Th.prime(h1)*E1
        g.W1 = -(t(x) %*% deltaE1)/length(E)
        if (!is.null(monotone))
            g.W1[monotone,] <- g.W1[monotone,]*W1[monotone,]
        gradient <- c(g.W1, g.W2)
    } else{
        y2 <- Th(h2)
        aug.y2 <- cbind(y2, 1)
        h3 <- aug.y2 %*% W3
        y3 <- To(h3)
        E <- y-y3
        deltaE3 <- 2*To.prime(h3)*E
        g.W3 <- -(t(aug.y2) %*% deltaE3)/length(E)
        if (!is.null(monotone))
            g.W3[1:(nrow(W3)-1),] <- g.W3[1:(nrow(W3)-1),]*W3[1:(nrow(W3)-1),]
        E2 <- deltaE3 %*% t(W3[1:(nrow(W3)-1),,drop=FALSE])
        deltaE2 <- Th.prime(h2)*E2
        g.W2 <- -(t(aug.y1) %*% deltaE2)/length(E)
        if (!is.null(monotone))
            g.W2[1:(nrow(W2)-1),] <- g.W2[1:(nrow(W2)-1),]*W2[1:(nrow(W2)-1),]
        E1 <- deltaE2 %*% t(W2[1:(nrow(W2)-1),,drop=FALSE])
        deltaE1 <- Th.prime(h1)*E1
        g.W1 <- -(t(x) %*% deltaE1)/length(E)
        if (!is.null(monotone))
            g.W1[monotone,] <- g.W1[monotone,]*W1[monotone,]
        gradient <- c(g.W1, g.W2, g.W3)
    }
    cost <- sum(E^2)/length(E)
    attr(cost, "gradient") <- gradient
    cost
}

