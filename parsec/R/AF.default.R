AF.default <-
function(y, z, w=rep(1, ncol(y)), k=sum(w), freq=rep(1, nrow(y)), ...) {
    res <- list(y = y, freq=freq)
    res$d <- sum(w)
    res$n <- sum(freq)
    res$z <- z
    res$w <- w
    res$k <- k
    res$rho <- function(x) x < z
    res$rho_k <- function(x) sum((x < z)*w) >= k
    res$g0 <- t(apply(y, 1, function(x) x < z))
    res$c <- apply(res$g0, 1, function(x) sum(x * w))
    res$Z_k <- res$c >= k
    res$q <- sum(res$Z_k * freq)
    res$H <- res$q/res$n
    res$A <- sum(res$c * res$Z_k * freq)/res$d/res$q
    res$M0 <- res$H*res$A
    class(res) <- "ophi"
    return(res)
}
