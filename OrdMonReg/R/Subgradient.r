Subgradient <- function(b, g1, w1, g2, w2, B, Gsi){

n <- length(g1)
grad <- rep(0, n)

# determine s_1, ..., s_k
sik <- matrix(NA, nrow = n, ncol = n)
for (i in 1:n){
    s <- 1:i
    if(i < n){tmp1 <- pmax(Gsi[i, s], b[s])}           
    if(i == n){tmp1 <- pmax(Gsi[i, s], c(b, B)[s])}    
    tmp2 <- s[min(tmp1) == tmp1]
    sik[i, 1:length(tmp2)] <- tmp2
} # end i

# compute subgradient
for (i in 1:n){ #change
    sis <- sik[i, ]
    sis <- sis[is.na(sis) == FALSE]
    Gisj <- Gsi[i, sis]
    if(i < n){gtype <- (Gisj >= b[sis])}  
    if((i == n) && (max(sis) == n)){gtype <- c((Gisj[-length(sis)] >= b[sis[-length(sis)]]), TRUE)}
    if((i == n) && (max(sis) < n)){gtype  <- (Gisj >= b[sis])}
    
    if(all(gtype) == FALSE){
            sim <- max(sis[gtype == FALSE])
            unitvec <- rep(0, n)
            unitvec[sim] <- 1
            grad <- grad + (2 * (max(Gsi[i, sim], b[sim]) - g1[i]) * w1[i]) * unitvec
    } # end if   
} # end i

# add second term to grad
grad <- grad + 2 * (c(b, B) - g2[1:n]) * w2[1:n] * c(rep(1, n-1), 0)

res <- list("grad" = grad)
return(res)
}

