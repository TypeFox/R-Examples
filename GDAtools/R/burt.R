burt <- function(data) {
    disj <- dichotom(data,out='numeric')
    res <- as.matrix(t(disj)) %*% as.matrix(disj)
    return(res)
    }
