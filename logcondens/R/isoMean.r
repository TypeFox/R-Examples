isoMean <- function(y, w){
    n <- length(y)
    k <- 1:n * 0
    gew <- 1:n * 0
    ghat <- 1:n * 0
    c <- 1
    k[c] <- 1
    gew[c] <- w[1]
    ghat[c] <- y[1]
    for (j in 2:n){
        c <- c + 1
        k[c] <- j
        gew[c] <- w[j]
        ghat[c] <- y[j]
        while (c >= 2 && ghat[max(1, c - 1)] >= ghat[c]){
            neu <- gew[c] + gew[c - 1]
            ghat[c - 1] <- ghat[c - 1] + (gew[c]/neu) * (ghat[c] - ghat[c - 1])
            gew[c - 1] <- neu
            c <- c - 1
        }
    }
    while (n >= 1){
        for (j in k[c]:n){ghat[j] <- ghat[c]}
        n <- k[c] - 1
        c <- c - 1
    }
    return(ghat)
}
