## -*- truncate-lines: t; -*-
## Time-stamp: <2014-03-26 17:12:19 CET (es)>

colSubset <- function(x) {

    nr <- dim(x)[1L]
    nc <- dim(x)[2L]    
    QR <- qr(x)
    
    if ((qrank <- QR$rank) == nc) {
        ans <- list(columns = seq_len(nc),
                    multiplier = diag(nc))
    } else {
        cols  <- QR$pivot[seq_len(qrank)]
        cols_ <- QR$pivot[(qrank + 1L):nc]
        D <- array(0, dim = c(qrank, nc))
        D[cbind(seq_along(cols), cols)] <- 1
        D[ ,cols_] <- qr.solve(x[ ,cols], x[ ,cols_])
        ans <- list(columns = cols,
                    multiplier = D)
    }
    ans
}
