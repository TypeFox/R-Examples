falk <- function(est, ks = NA){
#
# Calculate Falk's tail index estimator (aka negative Hill estimator), for
# an unknown endpoint.
#
# Input:
# - res : dlc object as output by logConDens in logcondens
#
# Kaspar Rufibach, 2006
#
n <- est$n
x <- est$xn        
        
# calculate quantiles 
c <- 1:n
q <- logcondens::quantilesLogConDens(ps = c / n, est)[, "quantile"]

v1 <- 1:n * NA
v2 <- v1
k0 <- 3:(n - 1)

if (identical(NA, ks)){ks <- k0}
ks <- ks[(ks %in% k0)]

for (k in ks){
    j <- 2:k
    
    # Falks based on quantiles of log-concave 
    temp <- (q[n] - q[n - j + 1]) / (q[n] - q[n - k])
    v1[k] <- sum(log(temp)) / (k - 1)
    
    # Falks based on order statistics 
    temp <- (x[n] - x[n - j + 1]) / (x[n] - x[n - k])
    v2[k] <- sum(log(temp)) / (k - 1)
}

res <- cbind("k" = 1:n, "logcon" = v1, "order" = v2) 

return(res)
}

