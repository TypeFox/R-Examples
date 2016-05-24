thinPNO <- function(x, n = 100){
	
	r <- range(x[, 1])
    r <- seq(from = r[1], to = r[2], length.out = n + 1)

    bd <- NULL
    for (i in 1:(length(r) - 1))
        bd <- rbind(bd, r[c(i, i + 1)])
    space <- bd[1, 2] - bd[1, 1]

    out <- matrix(nrow = n, ncol = dim(x)[2])
    out[, 1] <- bd[,1] + space/2
    colnames(out) <- colnames(x)
    for (i in seq(along = out[, 1])){
	    id <- which(x[, 1] > bd[i, 1] & x[, 1] <= bd[i, 2])
	    out[i, -1] <- apply(x[id, -1], 2, sum)
    }
    out
}