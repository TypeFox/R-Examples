cover2incidence.cover <-
function(g) {
    pow2 <- function(matrix) matrix %*% matrix
    g <- as.matrix(g)
    n <- dim(g)[1]
#     m <- floor(log(n-1, 2))  # per sicurezza potrei fare anche log(n, 2)
#     # ma floor dovrebbe essere sufficiente, pag: 210 Patil, Taillie
    m <- floor(log(n, 2))
    p <- g + diag(1, n)
    if(is.finite(m)) for(i in 1:m) p <- pow2(p)
    res <- p>0
    class(res) <- "incidence"
    return(res)
}
