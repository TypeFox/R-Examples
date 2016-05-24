`optim.fun6` <-
function (x, data, param,len,type) 
{
    L1 <- x[1]
    L2 <- x[2]
    L3 <- x[3]
    L4 <- x[4]
    data <- sort(data)
    q <- seq(0, 1, length = len)
    data.q <- quantile(data, prob = q,type=type)
    fitted.q <- qgl(q, L1, L2, L3, L4, param)
    not.inf<-!is.inf(fitted.q)

    check1 <- gl.check.lambda.alt(L1, L2, L3, L4, param)
    check2 <- fitted.q[1] <= min(data)
    check3 <- fitted.q[len] >= max(data)
    overall<-sum(c(check1,check2,check3),na.rm=T)
    
    if (overall==3) {
        # response <- sqrt(sum((data.q-fitted.q)^2))
        response <- sqrt(sum((data.q[not.inf]-fitted.q[not.inf])^2))
    }
    else if (overall!=3) {
        response <- NA
    }
    return(response)
}

