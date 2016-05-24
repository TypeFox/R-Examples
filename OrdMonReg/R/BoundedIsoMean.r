BoundedIsoMean <- function(y, w, a = NA, b = NA){
#
# Pool-adjacent-violaters-algorithm for a weighted isotonic
# mean, bounded by a lower bound a and an upper bound b.
#
# Input:
#   - y : data points
#   - w : corresponding weights
#   - a : vector of lower bound
#   - b : vector of upper bound
#
# Fadoua Balabdaoui and Kaspar Rufibach, October 2008

n <- length(y)
if (identical(a, NA)){a <- rep(-Inf, n)}
if (identical(b, NA)){b <- -a}
k <- 1:n * 0
ghat <- k
c <- 1
k[1] <- 1
ghat[1] <- MA(y, w, A = 1, a, b)    

for (j in 2:n){
    c <- c + 1
    k[c] <- j
    ghat[c] <- MA(y, w, A = j, a, b)

    while ((c >= 2) && (ghat[max(1, c - 1)] >= ghat[c])){

        ## use MA instead, general bounded case
        ind <- k[c - 1]:j
        ghat[c - 1] <- MA(y, w, A = ind, a, b)            
        c <- c-1
    } # end while
} # end for j

while (n >= 1){
    for (j in k[c]:n){ghat[j] <- ghat[c]}
    n <- k[c] - 1
    c <- c - 1
} # end while

return(ghat)
}

