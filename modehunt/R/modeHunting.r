modeHunting <- function(X.raw, lower = -Inf, upper = Inf, crit.vals, min.int = FALSE){

crit.vals <- unlist(crit.vals)

X <- preProcessX(X.raw, lower, upper)
n <- length(X)

sv <- sqrt(3 / (1 : (n - 2)))
cv <- sqrt(2 * (1 + log((n - 1) / (2 : n))))
Tv <- -cv
Dm <- NULL
Dp <- NULL
Dm.noadd <- NULL
Dp.noadd <- NULL

for (j in (n-2):1){
    tmp1 <- X[(j+1):n] - X[j]
    tmp2 <- cumsum(tmp1[1:(length(tmp1)-1)])
    tmp2 <- (2 * tmp2 / tmp1[2:length(tmp1)] - 1:(n-j-1)) * sv[1:(n-j-1)]
    
    ## Update of Tv:
    Tv[1:(n-j-1)] <- pmax(Tv[1:(n-j-1)], abs(tmp2) - cv[1:(n-j-1)])
    
    ## Update of Dm:
    JJ <- (-tmp2 - cv[1:(n-j-1)]) > crit.vals[1]
    JJ <- 1:length(JJ) * JJ
    JJ <- JJ[JJ>0]
    
    if (length(JJ) != 0){
        if (length(Dm[, 1]) == 0){
            Dm <- matrix(c(X[j], X[j + 1 + min(JJ)]), ncol = 2)} else {
            Dm <- rbind(c(X[j], X[j + 1 + min(JJ)]), Dm)}
    }
    
    ## Update of Dp:
    JJ <- (tmp2 - cv[1:(n-j-1)]) > crit.vals[1]
    JJ <- 1:length(JJ) * JJ
    JJ <- JJ[JJ>0]
    
    if (length(JJ) != 0){
        if (length(Dp[, 1]) == 0){
            Dp <- matrix(c(X[j], X[j + 1 + min(JJ)]), ncol =2)} else { 
            Dp <- rbind(c(X[j], X[j + 1 + min(JJ)]), Dp)}
    }
    
    ## Update of Dm.noadd:
    JJ <- -tmp2 > crit.vals[2]
    JJ <- 1:length(JJ) * JJ
    JJ <- JJ[JJ>0]
    
    if (length(JJ) != 0){
        if (length(Dm.noadd[, 1]) == 0){
            Dm.noadd <- matrix(c(X[j], X[j + 1 + min(JJ)]), ncol = 2)} else {
            Dm.noadd <- rbind(c(X[j], X[j + 1 + min(JJ)]), Dm.noadd)}
    }

    ## Update of Dp.noadd:
    JJ <- tmp2 > crit.vals[2]
    JJ <- 1:length(JJ) * JJ
    JJ <- JJ[JJ>0]
    
    if (length(JJ) != 0){
        if (length(Dp.noadd[, 1]) == 0){
            Dp.noadd <- matrix(c(X[j], X[j + 1 + min(JJ)]), ncol =2)} else { 
            Dp.noadd <- rbind(c(X[j], X[j + 1 + min(JJ)]), Dp.noadd)}
    }
}

## if min.int == TRUE, compute minimal intervals among all intervals
Dm.all <- NULL; Dp.all <- NULL; Dm.noadd.all <- NULL; Dp.noadd.all <- NULL

if (length(Dm[, 1]) > 0){
    Dm.all <- Dm[order(Dm[, 1]), ]
    if (min.int == TRUE){Dm.all <- minimalIntervals(Dm.all)}
} 

if (length(Dp[, 1]) > 0){
    Dp.all <- Dp[order(Dp[, 1]), ]
    if (min.int == TRUE){Dp.all <- minimalIntervals(Dp.all)}
} 

if (length(Dm.noadd[, 1]) > 0){
    Dm.noadd.all <- Dm.noadd[order(Dm.noadd[, 1]), ]
    if (min.int == TRUE){Dm.noadd.all <- minimalIntervals(Dm.noadd.all)}
} 

if (length(Dp.noadd[, 1]) > 0){
    Dp.noadd.all <- Dp.noadd[order(Dp.noadd[, 1]), ]
    if (min.int == TRUE){Dp.noadd.all <- minimalIntervals(Dp.noadd.all)}
} 

res <- list("Dp" = Dp.all, "Dm" = Dm.all, "Dp.noadd" = Dp.noadd.all, "Dm.noadd" = Dm.noadd.all)
return(res)
}


