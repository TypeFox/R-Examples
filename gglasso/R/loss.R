hubercls <- function(r, delta) {
    d <- matrix(NA, nrow(r), ncol(r))
    for (i in 1:nrow(r)) {
        for (j in 1:ncol(r)) {
            if (r[i, j] > 1) 
                d[i, j] <- 0 else {
                if (r[i, j] <= (1 - delta)) 
                  d[i, j] <- 1 - r[i, j] - 0.5 * delta else d[i, j] <- 0.5 * (1 - r[i, j])^2/delta
            }
        }
    }
    d
} 
