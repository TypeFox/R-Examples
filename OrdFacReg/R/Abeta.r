Abeta <- function(V, beta){
    q <- dim(V)[2]
    A <- (1:q)[t(V) %*% beta >= 0]
    return(A)
    }
