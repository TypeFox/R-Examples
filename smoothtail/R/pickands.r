pickands <- function(est, ks = NA){
#
# Calculate Pickands' tail index estimator, for an unknown endpoint.
#
# Input:
# - x : Vector of observations
#
# Kaspar Rufibach, 2010
#
n <- est$n
x <- est$xn 
res <- NA
        
# calculate quantiles 
c <- 1:n
q <- logcondens::quantilesLogConDens(ps = c / n, est)[, "quantile"]

v1 <- 1:n * NA
v2 <- v1 

k0 <- 4:n
if (identical(NA, ks)){ks <- k0}
ks <- ks[(ks %in% k0)]

for (k in ks){
    k2 <- floor(k / 4)
    
    # Pickands based on quantiles of log-concave 
    q1 <- logcondens::quantilesLogConDens((n - k / 4 + 1) / n, est)[, "quantile"]
    q2 <- logcondens::quantilesLogConDens((n - k / 2 + 1) / n, est)[, "quantile"]
    q3 <- logcondens::quantilesLogConDens((n - k     + 1) / n, est)[, "quantile"]
    v1[k] <- (q1 - q2) / (q2 - q3)

    # Pickands based on order statistics 
    v2[k] <- (x[n - k2 + 1] - x[n - 2 * k2 + 1]) / (x[n - 2 * k2 + 1] - x[n - 4 * k2 + 1])
}

v1 <- log(v1) / log(2)
v2 <- log(v2) / log(2)

# output result
res <- cbind("k" = 1:n, "logcon" = v1, "order" = v2) 
return(res)
}

