"compdelta" <-
function(Pi){
    #   From Zucchini Sydney notes
    #   Calculate delta assuming stationarity
    m <- ncol(Pi)
    a <- t(Pi) - diag(1, m, m)
    a[m,] <- 1
    b <- c(rep(0, (m-1)), 1)
    delta <- solve(a, b)
    return(delta)
}

