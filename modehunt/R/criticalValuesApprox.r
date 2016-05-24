criticalValuesApprox <- function(n, d0 = 2, m0 = 10, fm = 2, alpha = 0.05, gam = 2, tail = 10, M = 10 ^ 5, display = 0, path = NA){

fd <- sqrt(fm)
sv <- sqrt(1:n / 3)
del <- (1:(n + 1)) / (n + 2) 
cv <- sqrt(2 * (1 - log(del)))

n.blocks <- floor(log(n / m0) / log(fm))

Ts.mat.add <- rep(0, M) - sqrt(2)
Ts.mat.noadd <- Ts.mat.add
Ts.mat.block <- matrix(0, nrow = M, ncol = n.blocks) - sqrt(2)

for (sim in 1:M){
    
    U <- sort(runif(n, min = 0, max = 1))
    cU <- cumsum(U)
    d <- d0
    m <- m0
    block <- n.blocks
    
    for (block in 1:n.blocks){

        d <- myRound(d0 * fd ^ (n.blocks - block)) 
        m <- myRound(m0 * fm ^ (n.blocks - block))
        
        minspan <- d * ceiling((m + 1) / d)   
        maxspan <- d * floor(fm * m / d)
        
        for (j in seq(1, n - m + 1, by = d)){
            k <- seq(min(n, j + minspan), min(n, (j + maxspan)), by = d)
            if (length(k) > 0){
                Tjk <- (2 / (U[k] - U[j]) * (cU[k - 1] - cU[j] - U[j] * (k - j - 1)) - (k - j - 1))        
                Tjk <- abs(Tjk) / sv[k - j - 1]
                Ts.mat.add[sim] <- max(Ts.mat.add[sim], max(Tjk - cv[k - j]))
                Ts.mat.noadd[sim] <- max(Ts.mat.noadd[sim], max(Tjk))        
                Ts.mat.block[sim, block] <- max(Ts.mat.block[sim, block], max(Tjk))
            } #if
        } #j, k
    } #for
    
    if ((sim / 100 == myRound(sim / 100)) & (display == 1)){
        print(paste("n = ", n, " / sim = ", sim, sep = ""))
        name <- paste(path, "critvals_approx_alpha=", alpha, "_n=", n, ".txt", sep = "")
        if (is.na(path) == 0){write.table(sim, file = name, row.names = FALSE, col.names = FALSE)} 
    }             
} #for

# compute critical values for each block
# raw search for \tilde alpha
Ts.mat.block.sort <- apply(Ts.mat.block, 2, sort)
cutoffs <- 1 : n.blocks 

# now start bisection to find right tail:
i.l <- 1 
i.r <- M
while (i.r - i.l > 1){
  i <- myRound((i.r + i.l) / 2)
  for (c in 1 : n.blocks){cutoffs[c] <- Ts.mat.block.sort[M - myRound((M - i) * (1 + tail) ^ gam / ((c + tail) ^ gam)), c]}
  count <- 0 
  for (d in 1 : M){count <- count + (sum(Ts.mat.block[d, ] >= cutoffs) > 0)}
  if ((count / M) > alpha){i.l <- i} else {i.r <- i}
}

# return critical values
quan <- myRound((1 - alpha) * M)
Ts.mat.add <- sort(Ts.mat.add)
Ts.mat.noadd <- sort(Ts.mat.noadd)

withadd <- Ts.mat.add[quan]
noadd <- Ts.mat.noadd[quan]

# generate output
return(list("approx" = c(withadd, noadd), "block" = cutoffs))
}




