getBounds <-
function(x, dec=2) {
    if(!is.numeric(x))   { stop("x must be numeric") }
    if(!is.numeric(dec)) { stop("dec must be numeric") }
    if(dec %% 1)         { stop("dec must be an integer") }

    adj <- 10^dec
    return(c(floor(min(x)*adj), ceiling(max(x)*adj)) / adj)
}

getMaxNorm <-
function(x, dec=2) {
    if(!is.numeric(x)) { stop("x must be numeric") }

    mm <- range(x)                       # minimum and maximum
    return(getBounds(dnorm(seq(mm[1], mm[2], length.out=200),
                           mean(x), sd(x)), dec))
}
