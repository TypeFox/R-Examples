"coefZ" <-
function(X){
    X <- check.data(X)
    N <- nrow(X)
    S <- var(X)
    Sij <- outer(apply(X,2,var),apply(X,2,var),"*")
    Zij <- (S * sqrt(N-1))/sqrt(Sij)
    diag(S) <- diag(Sij) <- diag(Zij) <- 0
    Zi <- (apply(S,1,sum) * sqrt(N-1))/ sqrt(apply(Sij,1,sum))       
    Z  <- (sum(S)/2 * sqrt(N-1))/ sqrt(sum(Sij)/2)
    return(list(Zij=Zij,Zi=Zi,Z=Z))
}

