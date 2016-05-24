"coefH" <-
function(X){
    X <- check.data(X)
    S <- var(X)
    Smax <- var(apply(X, 2, sort))
    Hij <- S/Smax
    diag(S) <- 0
    diag(Smax) <- 0
    Hi <- apply(S, 1, sum)/apply(Smax, 1, sum)
    H <- sum(S)/sum(Smax)
    return(list(Hij=Hij,Hi=Hi,H=H))
}

