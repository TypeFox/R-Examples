generalizedPick <- function(est, c = 0.75, gam0 = -1 / 3, ks = NA){

## estimators developped in Segers (2005), with the mixing measure
## given in Theorem 4.1, (i)
## note that the j in Theorem 4.1 does NOT correspond to the j in 
## Definition 2.1 of the paper.
n <- est$n
x <- est$xn 

delta <- abs(gam0 + 1 / 2) - 1 / 2
res <- NA
w <- 1:n
q <- logcondens::quantilesLogConDens(ps = w / n, est)[, "quantile"]
v1 <- 1:n * NA
v2 <- v1

k0 <- 1 : (n - 1)
if (identical(NA, ks)){ks <- k0}
ks <- ks[(ks %in% k0)]

for (k in ks){

    # smoothed observations
    j <- 1:k  
    lambda1 <- lambdaGenPick(t = j / k, delta, c)
    lambda2 <- lambdaGenPick(t = (j - 1) / k, delta, c)        
    v1[k] <- sum((lambda1 - lambda2) * log(q[n - floor(c * j)] - q[n - j]))
    
    # original observations
    v2[k] <- sum((lambda1 - lambda2) * log(x[n - floor(c * j)] - x[n - j]))
    }

res <- cbind(k = 1:n, logcon = v1, order = v2)    
return(res)
}
