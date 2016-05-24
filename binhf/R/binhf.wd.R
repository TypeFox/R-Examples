binhf.wd <-
function (x, binsize = 1, print.info = FALSE) 
{
    den <- NULL
    f <- NULL
    xwd <- ht(x)
    c <- xwd$s
    d <- xwd$d
    n <- length(d) + 1
    cnew <- c
    dnew <- d
    J <- log2(n)
    number <- binsize
    den <- cnew * (number - cnew)/(number)
    f <- dnew
    r <- (den > 0)
    q <- !r
    f[r] <- dnew[r]/sqrt(den[r])
    f[q] <- 0
    transformed <- ht.inv(c(d, cnew[length(cnew)])) 
      
    if (print.info == TRUE) {
        cat("cnew:\n")
        print(cnew)
        cat("d:\n")
        print(d)
        cat("den:\n")
        print(den)
        cat("f:\n")
        print(f)
    }
    if (is.list(transformed)) {
        transformed <- transformed$res
    }
    l <- list(transformed = transformed, cnew = cnew)
    return(l)
}

