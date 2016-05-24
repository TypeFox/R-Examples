fix.theta <-
function (theta, Q) 
{
    M <- t(theta) %*% Q %*% theta
    eM <- eigen(M, symmetric = TRUE)
    scale(theta %*% eM$vectors, FALSE, sqrt(eM$values))
}

