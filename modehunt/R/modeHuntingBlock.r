modeHuntingBlock <- function(X.raw, lower = -Inf, upper = Inf, d0 = 2, m0 = 10, fm = 2, crit.vals, min.int = FALSE){

crit.vals <- unlist(crit.vals)

X <- preProcessX(X.raw, lower, upper)
n <- length(X)
cX <- cumsum(X)
sv <- sqrt((1 : (n - 2)) / 3)
fd <- sqrt(fm)

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

            Tjks <- rbind(Tjks, cbind(j * rep(1, length(k)), k, Tjk, 1 / sv[k - j - 1]))        
        } ## if
    } ## end j                  
} ## end block in 1:n.blocks
    
Tjks <- Tjks[order(Tjks[, 1]), ]
    
## add number of observations as sixth column and sort accordingly
Tjks <- cbind(Tjks, Tjks[, 2] - Tjks[, 1] + 1)
Tjks <- Tjks[order(Tjks[, 5]), ]

## block sizes
bs <- blocks(n, m0, fm)
Dp.all <- NULL; Dm.all <- NULL

for (block in 1:n.blocks){
    ints <- Tjks[Tjks[, 5] >= bs[block, 1], ]
    ints <- ints[ints[, 5] <= bs[block, 2], ]
    n.ints <- length(ints[, 1])

    ## rescale critical values
    cjk <- ints[, 4] ^ (-1) * crit.vals[block]
    
    ## find intervals with significant increase within current block
    Dp <- ints[ints[, 3] >= cjk, 1:2]

    ## find intervals with significant decrease within current block
    Dm <- ints[-ints[, 3] >= cjk, 1:2]
       
    Dp.all <- rbind(Dp.all, Dp)
    Dm.all <- rbind(Dm.all, Dm)  
} ## for

## if min.int == TRUE, compute minimal intervals among all intervals
if (length(Dp.all[, 1]) > 0){
    Dp.all <- Dp.all[order(Dp.all[, 1]), ]
    if (min.int == TRUE){Dp.all <- minimalIntervals(Dp.all)}
    Dp.all <- cbind(X[Dp.all[, 1]], X[Dp.all[, 2]])
} else {Dp.all <- NULL}

if (length(Dm.all[, 1]) > 0){
    Dm.all <- Dm.all[order(Dm.all[, 1]), ]
    if (min.int == TRUE){Dm.all <- minimalIntervals(Dm.all)}
    Dm.all <- cbind(X[Dm.all[, 1]], X[Dm.all[, 2]])
} else {Dm.all <- NULL}

## generate output
res <- list("Dp" = Dp.all, "Dm" = Dm.all)
return(res)
}
