binaryRep <-
function(data, m = .Machine$double.digits) {
    x <- data
    n <- length(x)
    
    xSign <- sign(x)
    x <- xSign * x
    exponent <- ifelse(x > 0, floor(1+log(x, 2)), 0)
    x <- x/2 ^ exponent
    pwrs <- binaryRepPowers(n, m)
    x <- matrix(x, n, m)
    xpwrs <- x * pwrs
    xrep <- trunc(xpwrs)/pwrs
    bits <- (xrep %*% binaryRepA(m)) *pwrs
    
    bits[] <- as.integer(bits[])
    new("binaryRep", original = data,
        sign = as.integer(xSign),
        exponent = as.integer(exponent),
        bits = binaryRepBits(bits))
}
