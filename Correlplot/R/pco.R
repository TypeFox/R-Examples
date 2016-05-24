`pco` <-
function (Dis) {  
    n <- nrow(Dis)
    p <- ncol(Dis)
    A <- -0.5 * (Dis * Dis)
    n <- nrow(Dis)
    I <- diag(rep(1, n))
    H <- I - (1/n) * rep(1,n)%*%t(rep(1,n))
    B <- H %*% A %*% H
    Out <- eigen(B)
    Dl <- diag(Out$values)
    k <- sum(Out$values >= 0)
    cat("There are ", k, " non-negative eigenvalues.\n")
    Dk <- Dl[1:k, 1:k]
    V <- Out$vectors
    Vk <- V[, 1:k]
    if (k == 1) 
        PC <- Vk * Dk
    else PC <- Vk %*% sqrt(Dk)
    Total <- sum(Dl)
    ev <- diag(Dl)
    fr <- diag(Dl)/Total
    cu <- cumsum(fr)
    decom <- cbind(ev, fr, cu)
    print(round(decom, digits = 2))
    return(list(PC = PC, Dl = Dl, Dk = Dk, B = B, decom = decom))
}
