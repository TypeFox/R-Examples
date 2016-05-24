modeHuntingApprox <- function(X.raw, lower = -Inf, upper = Inf, d0 = 2, m0 = 10, fm = 2, crit.vals, min.int = FALSE){

crit.vals <- unlist(crit.vals)

X <- preProcessX(X.raw, lower, upper)
n <- length(X)
cX <- cumsum(X)
fd <- sqrt(fm)

# compute test statistics
sv <- sqrt((1 : (n - 2)) / 3)
del <- (1: (n + 1)) / (n + 2) 
cv <- sqrt(2 * (1 - log(del)))

n.blocks <- floor(log(n / m0) / log(fm))
Tjks <- NULL 

for (block in 1:n.blocks){

    d <- myRound(d0 * fd ^ (n.blocks - block)) 
    m <- myRound(m0 * fm ^ (n.blocks - block))
    
    minspan <- d * ceiling((m + 1) / d)   
    maxspan <- d * floor(fm * m / d)
  
    for (j in seq(1, n - m + 1, by = d)){
        k <- seq(min(n, j + minspan), min(n, (j + maxspan)), by = d)
        if (length(k) > 0){
            Tjk <- (2 / (X[k] - X[j]) * (cX[k - 1] - cX[j] - X[j] * (k - j - 1)) - (k - j - 1))

            Tjks <- rbind(Tjks, cbind(j * rep(1, length(k)), k, Tjk, 1 / sv[k - j - 1], cv[k - j]))         
        } ## if
    } ## end j                  
} ## end block in 1:n.blocks
    
Tjks <- Tjks[order(Tjks[, 1]), ]

# choose decisions on H1, with additive correction (use crit.vals(1) as critical value)
cjk <- Tjks[, 4] ^ (-1) * (Tjks[, 5] + crit.vals[1])
Dp <- Tjks[Tjks[, 3] >= cjk, 1:2]
Dm <- Tjks[-Tjks[, 3] >= cjk, 1:2]

# choose decisions on H1, without additive correction (use crit.vals(2) as critical value)
cjk <- Tjks[, 4] ^ (-1) * crit.vals[2]
Dp.noadd <- Tjks[Tjks[, 3] >= cjk, 1:2]
Dm.noadd <- Tjks[-Tjks[, 3] >= cjk, 1:2]

# if min.int == TRUE, compute minimal intervals among all intervals
Dm.all <- NULL; Dp.all <- NULL; Dm.noadd.all <- NULL; Dp.noadd.all <- NULL

if (length(Dp[, 1]) > 0){
    Dp <- Dp[order(Dp[, 1]), ]
    if (min.int == TRUE){Dp <- minimalIntervals(Dp)}
    Dp.all <- cbind(X[Dp[, 1]], X[Dp[, 2]])
} 
    
if (length(Dm[, 1]) > 0){
    Dm <- Dm[order(Dm[, 1]), ]
    if (min.int == TRUE){Dm <- minimalIntervals(Dm)}
    Dm.all <- cbind(X[Dm[, 1]], X[Dm[, 2]])
}

if (length(Dp.noadd[, 1]) > 0){
    Dp.noadd <- Dp.noadd[order(Dp.noadd[, 1]), ]
    if (min.int == TRUE){Dp.noadd <- minimalIntervals(Dp.noadd)}
    Dp.noadd.all <- cbind(X[Dp.noadd[, 1]], X[Dp.noadd[, 2]])
}
    
if (length(Dm.noadd[, 1]) > 0){
    Dm.noadd <- Dm.noadd[order(Dm.noadd[, 1]), ]
    if (min.int == TRUE){Dm.noadd <- minimalIntervals(Dm.noadd)}
    Dm.noadd.all <- cbind(X[Dm.noadd[, 1]], X[Dm.noadd[, 2]])
}

# generate output
res <- list("Dp" = Dp.all, "Dm" = Dm.all, "Dp.noadd" = Dp.noadd.all, "Dm.noadd" = Dm.noadd.all)
return(res)
}

