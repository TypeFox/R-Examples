gibbs_admkr_erro = function(xh, inicost, k, errorsizp, errorprob, data_x, data_y)
{
    dm = ncol(data_x)
    fx = inicost
    dv = rnorm(1) * errorsizp
    xh[dm + 1] = xh[dm + 1] + dv
    fy = cost_admkr(xh, data_x, data_y)
    r = fx - fy
    if(class(r) != "numeric")
    {
        stop("Please re-run the code.")
    }
    if(r > 0)
    {
        accept = 1
    }
    else
    {
        if(runif(1) < exp(r))
        {
            accept = 1
        }
        else
        {
            accept = 0
        }
    }
    c = errorsizp/(errorprob * (1 - errorprob))
    if(accept == 1)
    {
        accept_erro = 1
        inicost = fy
        errorsizp = errorsizp + c * (1 - errorprob)/k
    }
    else
    {
        accept_erro = 0
        xh[dm+1] = xh[dm+1] - dv
        errorsizp = errorsizp - c * errorprob/k
    }
    return(list(x = xh, cost = inicost, accept_erro = accept_erro, errorsizp = errorsizp))
}

