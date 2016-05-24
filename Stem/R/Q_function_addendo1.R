`Q_function_addendo1` <-
function(sigma2omega, n, Sigmastar, B){
        Q_addendo1  <- n*log(det(sigma2omega * Sigmastar)) + sum(diag(solve(sigma2omega * Sigmastar)%*% B))
        return(Q_addendo1)
}

