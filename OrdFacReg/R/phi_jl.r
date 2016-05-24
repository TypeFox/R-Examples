phi_jl <- function(p, f, kj){

# compute inequality nr for given j, l

c <- p - f
kj0 <- c(rep(0, c), kj)
mjls <- NULL

for (j in (c+1):p){
    for (l in 2:(kj0[j]-1)){
        h <- (c+1):j
        i <- sum(kj0[h-1])+(l-1) - 2 * (j-c-1)
        mjls <- rbind(mjls, c(j, l, i))}
}
dimnames(mjls) <- list(NULL, c("j", "l", "i"))
return(mjls)
}


## neue Version (fuer >= 0 setting)
phi_jl <- function(p, f, kj){

# compute inequality nr for given j, l

c <- p - f
kj0 <- c(rep(0, c), kj)
mjls <- NULL

for (j in (c+1):p){
    for (l in 2:(kj0[j])){
        h <- (c+1):j
        i <- sum(kj0[h - 1]) + (l - 1) - (j - c - 1)
        mjls <- rbind(mjls, c(j, l, i))}
}
dimnames(mjls) <- list(NULL, c("j", "l", "i"))
return(mjls)
}
