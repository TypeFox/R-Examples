falkMVUE <- function(est, omega, ks = NA){
#
# Calculate Falk's minimum variance unbiased tail index estimator, for a known endpoint.
#
# Input:
# - x   : Vector of quantiles. If the ordinary Pickand estimator is to be calculated,
#         x simply equals the vector of ordered observations.
# - gam : Known tail index.
#
# Kaspar Rufibach, 2010
#
n <- est$n
x <- est$xn 
v1 <- 1:n*NA
v2 <- v1

if (omega < x[n]){
    cat("omega must be greater than max(x)!")} else {
        
    # calculate quantiles 
    c <- 1:n
    q <- logcondens::quantilesLogConDens(ps = c / n, est)[, "quantile"]

    k0 <- 2:(n - 1)
    if (identical(NA, ks)){ks <- k0}
    ks <- ks[(ks %in% k0)]

    for (k in ks){
        j <- 2:k
    
        # Falks based on quantiles of log-concave 
        temp <- (omega - q[n-j+1]) / (omega - q[n-k])
        v1[k] <- sum(log(temp)) / k
    
        # Falks based on order statistics 
        temp <- (omega - x[n-j+1]) / (omega - x[n-k])
        v2[k] <- sum(log(temp)) / k
    }

res <- cbind("k" = 1:n, "logcon" = v1, "order" = v2) 

return(res)}
}

