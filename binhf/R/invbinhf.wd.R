invbinhf.wd <-
function (transformed, binsize = 1, print.info = FALSE) 
{
    corig <- transformed$cnew
    transformed <- transformed$transformed
    if (is.list(transformed)) {
        transformed <- transformed$res
    }
    c <- ht2(transformed)$s
    d <- ht2(transformed)$d
    n <- length(d) + 1
    den <- f <- NULL
    dnew <- d
    cnew <- c
    number <- binsize
    den <- (cnew * (number - cnew)/(number))
    f <- dnew
    r <- (den > 0)
    q <- !r
    f[r] <- dnew[r] * sqrt(den[r])
    f[q] <- 0
        estimate <- hf.inv2(list(d = d, s = cnew), binsize = binsize) 
    if (print.info == TRUE) {
        cat("c is\n")
        print(c)
        cat("cnew is\n")
        print(cnew)
        cat("dnew is\n")
        print(dnew)
        cat("den is \n")
        print(den)
        cat("f is\n")
        print(f)
    }
    if (is.list(estimate)) {
        estimate <- estimate$res
    }
    return(estimate)
}

