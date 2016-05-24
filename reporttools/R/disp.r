disp <- function(n, d1 = 2, d2 = 1){
    n <- as.numeric(n)
    ind.na <- is.na(n) == FALSE
    ind <- (abs(n) >= 10^-d1)
    n[ind.na & ind] <- format(round(as.numeric(n[ind.na & ind]), d1), nsmall = d1)
    tmp <- n[(ind.na & ind) == FALSE]
    for (i in 1:length(tmp)){tmp[i] <- format(as.numeric(tmp[i]), digits = d2, scientific = FALSE)}
    n[(ind.na & ind) == FALSE] <- tmp
    return(n)
}
