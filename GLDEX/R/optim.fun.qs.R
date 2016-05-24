`optim.fun.qs` <-
function (x, data, param1, param2,len,type) 
{
    L1 <- x[1]
    L2 <- x[2]
    L3 <- x[3]
    L4 <- x[4]
    M1 <- x[5]
    M2 <- x[6]
    M3 <- x[7]
    M4 <- x[8]
    p <- x[9]
    check1 <- gl.check.lambda.alt(L1, L2, L3, L4, param1)
    check2 <- gl.check.lambda.alt(M1, M2, M3, M4, param2)
    if (check1 == TRUE && check2 == TRUE) {
        check3 <- min(qgl(0, L1, L2, L3, L4, param1), qgl(0, 
            M1, M2, M3, M4, param2)) <= min(data)
        check4 <- max(qgl(1, L1, L2, L3, L4, param1), qgl(1, 
            M1, M2, M3, M4, param2)) >= max(data)
    }
    else {
        check3 <- FALSE
        check4 <- FALSE
    }
    check5 <- p > 0
    check6 <- p < 1
    if(as.logical(check5 * check6)) {
        first.q <- seq(0, 1, length = len * p)
        second.q <- seq(0, 1, length = len * (1 - p))
        theo.quantile <- sort(c(qgl(first.q, L1, L2, L3, L4, param1), 
        qgl(second.q, M1, M2, M3, M4, param2)))
        fitted.q <- (theo.quantile[!is.inf(theo.quantile)])
}
if(as.logical(check1 * check2 * check3 * check4 * check5 * check6)) {
response <-  sqrt(sum((quantile(data,sort(c(first.q, 
second.q))[!is.inf(theo.quantile)],type=type)-fitted.q)^2))
}
else if(as.logical(check1 * check2 * check3 * check4 * check5 * check6) == 0) {
response <- NA
}
return(response)
}

