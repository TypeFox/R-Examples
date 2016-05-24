fdr.est <-
function(p)
{
    m <- length(ind <- which(!is.na(p)))
    fdr <- rep(NA, length(p))
    stat <- cbind(1:length(p), p, fdr)
    stat[ind, 3] <- unlist(lapply(stat[ind, 2], function(x) {
        c <- length(which(stat[ind, 2] <= x))
        m * x/c
    }))
    stat[ind, ] <- stat[ind, ][order(stat[, 2], decreasing = TRUE), 
        ]
    stat[ind, 3] <- cummin(stat[ind, 3])
    fdr <- stat[order(stat[, 1]), 3]
    fdr[which(fdr > 1)] <- 1
    return(fdr)
}
