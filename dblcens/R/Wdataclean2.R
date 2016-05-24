####################################################
## this is the same function as in emplik package ##
####################################################
Wdataclean2 <- function (z, d, wt = rep(1, length(z))) 
{
    niceorder <- order(z, -d)
    sortedz <- z[niceorder]
    sortedd <- d[niceorder]
    sortedw <- wt[niceorder]
    n <- length(sortedd)
    y1 <- sortedz[-1] != sortedz[-n]
    y2 <- sortedd[-1] != sortedd[-n]
    y <- y1 | y2
    ind <- c(which(y | is.na(y)), n)
    csumw <- cumsum(sortedw)
    list(value = sortedz[ind], dd = sortedd[ind], weight = diff(c(0, 
        csumw[ind])))
}

