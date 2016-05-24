lambdaGenPick <- function(t, delta, c = 0.75){

## function lambda_{delta, rho}^c given in Theorem 4.1
## of Segers (2005)

j <- ceiling(log(t) / log(c))  
res <- (1 - c ^ (1 + delta)) * t
if (delta != 0){res <- res * (1 - c ^ (delta * j)) / (1 - c ^ delta)}

res[t == 0] <- 0 
return(res)
}

