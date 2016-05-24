S.piPS <- function (n, x, e = runif(length(x))) 
{
N <- length(x)
x1 <- sort(x, decreasing = TRUE)    
Pik <- PikPPS(n, x1)    
    
    V <- cumsum(Pik)
    nk <- matrix(0, N, 1)
    d <- matrix(0, N, 1)
    I <- matrix(0, N, 1)
    sam <- matrix(0, N, 1)
    if (e[1] < Pik[1]) {
        I[1] <- 1
        sam[1] <- 1
    }
    for (k in 2:N) {
        nk[k] <- nk[k - 1] + I[k - 1]
        d[k] <- Pik[k] * (n - nk[k])/(n - V[k - 1])
        if (e[k] <= d[k]) {
            I[k] <- 1
            sam[k] <- cumsum(I[1:(k - 1)])[(k - 1)] + I[k]
        }
    }
samp <- rev(order(x))[which(sam != 0)]
Pik1 <- PikPPS(n, x) 
Pik.s <- Pik1[samp]
return(cbind(samp, Pik.s))
}
