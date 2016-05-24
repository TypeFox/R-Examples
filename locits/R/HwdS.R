HwdS <-
function (x) 
{
    n <- length(x)
    J <- IsPowerOfTwo(n)
    if (is.na(J)) 
        stop("Length of x is not a power of two")
    ans <- vector("list", J)
    for (j in (J - 1):0) {
        cf <- 2^(-(J - j)/2)
        f <- c(rep(-cf, 2^(J - j - 1)), rep(cf, 2^(J - j - 1)))
        ans[[j + 1]] <- getridofendNA(filter(x = x, filter = f, sides = 2))
    }
    ans[[1]] <- rep(0, n)
    answdS <- wd(rep(0, n), filter.number = 1, family = "DaubExPhase", 
        type = "station")
    for (j in (J - 1):0) answdS <- putD(answdS, lev = j, ans[[j + 
        1]])
    answdS
}
