alts2dups <-
function (x, ifalt = FALSE) 
{
    ndup <- length(x) / 2
    x1 <- numeric(ndup)
    x2 <- numeric(ndup)
    if (ifalt) {
        for (i in 1:ndup) {
            j <- 2 * (i - 1) + 1
            x1[i] <- x[j]
            x2[i] <- x[j + 1]
        }
    }
    else {
        for (i in 1:ndup) {
            x1[i] <- x[i]
            x2[i] <- x[ndup + i]
        }
    }
    xname<- deparse(substitute(x))
    x1name <- paste(xname, ".1", sep="")
    x2name <- paste(xname, ".2", sep="")
    xx <- matrix(c(x1, x2), nrow = ndup, ncol = 2,
                 dimnames = list(seq(1, ndup), c(x1name, x2name)))
    return(xx = xx)
}
