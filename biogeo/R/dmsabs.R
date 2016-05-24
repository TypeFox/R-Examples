dmsabs <-
function (coordstr, fmt) 
{
    x <- sub(" +$", "", coordstr)
    f1 <- substr(x, 1, 1)
    x <- ifelse(f1 == " ", sub(" ", "", x), x)
    f1 <- substr(x, 1, 1)
    x <- ifelse(f1 == " ", sub(" ", "", x), x)
    let <- c("WESN")
    L <- rep(0, 4)
    L[1] <- length(grep("W", toupper(x)))
    L[2] <- length(grep("E", toupper(x)))
    L[3] <- length(grep("S", toupper(x)))
    L[4] <- length(grep("N", toupper(x)))
    if (sum(L) == 0) {
        er <- "missing letters"
        letr <- "-"
        exclude <- 1
    }
    else {
        er <- "-"
        exclude <- 0
        f <- which(L == 1)
        letr <- substr(let, f, f)
        x <- sub(letr, "", x)
    }
    x <- sub(" +$", "", x)
    f1 <- substr(x, 1, 1)
    x <- ifelse(f1 == " ", sub(" ", "", x), x)
    g <- data.frame(coord = x, L = letr, error = er, exclude, 
        stringsAsFactors = F)
    str <- g$coord
    fdd <- unlist(gregexpr("d", fmt))
    dd <- as.numeric(substr(str, fdd[1], fdd[length(fdd)]))
    fmm <- unlist(gregexpr("m", fmt))
    mm <- as.numeric(substr(str, fmm[1], fmm[length(fmm)]))
    fss <- unlist(gregexpr("s", fmt))
    ss <- as.numeric(substr(str, fss[1], fss[length(fss)]))
    let <- ifelse(g$L == "-", NA, g$L)
    dms <- data.frame(dd, mm, ss, L = let)
    return(dms)
}
