post.prob <-
function (y, x, B, sigma, w1, lambda) {
     n <- nrow(y)
     maxold  <- matrix (-999999999,n,1)
     dif <- y - x%*%B
    wp1 <- matrix(log(w1),n,1)
    wp2 <- matrix(log(1-w1),n,1)
    s1 <- solve(sigma)
    s2 <- s1 / (lambda + 1)
    q1 <- matrix(tensorizza (dif, s1),n,1)
    q2 <- matrix(tensorizza (dif, s2),n,1)
    rm (s1,s2)
    gc()
    lntau1 <- -0.5*q1 + wp1
    lntau2 <- -0.5*q2 + wp2
    maxln <- pmax (maxold, lntau1, lntau2)
    maxold <- maxln
    appo1 <- lntau1 - maxold
    appo2 <- lntau2 - maxold
    sigma2 <- (1+lambda)* sigma
    
    dd1 <- sqrt(det(as.matrix(sigma)))
    dd2 <- sqrt(det(as.matrix(sigma2)))
    numtau1 <- exp(appo1)/dd1
    numtau2 <- exp(appo2)/dd2
    dentau <- (numtau1+numtau2)
    tau1 <-  numtau1 / dentau
    tau1
}

