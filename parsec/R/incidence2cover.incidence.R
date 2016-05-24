incidence2cover.incidence <-
function(z) {
    n <- dim(z)[1]
    eta <- z - diag(1, n)
    res <- (eta - 1*((eta%*%eta)>0))>0
    class(res) <- "cover"
    return(res)
}
