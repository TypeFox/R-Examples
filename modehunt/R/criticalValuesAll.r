criticalValuesAll <- function(n, alpha = 0.05, M = 10 ^ 5, display = 0, path = NA){

sv <- sqrt(1:n / 3)
del <- (1:(n + 1)) / (n + 2) 
cv <- sqrt(2 * (1 - log(del)))

Ts.mat.add <- rep(0, M) - sqrt(2)
Ts.mat.noadd <- Ts.mat.add

for (sim in 1:M){
    
    U <- sort(runif(n, min = 0, max = 1))
    cU <- cumsum(U)
    for (j in 1 : (n - 2)){
        k <- (j + 2) : n
        Tjk <- abs((2 / (U[k] - U[j]) * (cU[k - 1] - cU[j] - U[j] * (k - j - 1)) - (k - j - 1)) / sv[k - j - 1])
            
        Ts.mat.add[sim] <- max(Ts.mat.add[sim], max(Tjk - cv[k - j]))
        Ts.mat.noadd[sim] <- max(Ts.mat.noadd[sim], max(Tjk))        
    } #j

    if ((sim / 100 == myRound(sim / 100)) & (display == 1)){
        print(paste("n = ", n, " / sim = ", sim, sep = ""))
        name <- paste(path, "critvals_all_alpha=", alpha, "_n=", n, ".txt", sep = "")
        if (is.na(path) == 0){write.table(sim, file = name, row.names = FALSE, col.names = FALSE)}
    }
            
} #for

# return critical values
quan <- myRound((1 - alpha) * M)
Ts.mat.add <- sort(Ts.mat.add)
Ts.mat.noadd <- sort(Ts.mat.noadd)

withadd <- Ts.mat.add[quan]
noadd <- Ts.mat.noadd[quan]

## generate output
res <- c("withadd" = withadd, "noadd" = noadd)
return(res)
}

