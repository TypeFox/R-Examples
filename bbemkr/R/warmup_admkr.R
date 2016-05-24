warmup_admkr = function(x, inicost, mutsizp, errorsizp, warm = 100, prob = 0.234, errorprob = 0.44, data_x, data_y)
{
    mutsizpwarm = errorsizpwarm = vector(,warm)
    for(k in 1:warm)
    {
        dummy = gibbs_admkr_nw(xh = x, inicost = inicost, k = k, mutsizp = mutsizp, prob = prob, data_x = data_x, data_y = data_y)
        x = dummy$x
        inicost = dummy$cost
        mutsizpwarm[k] = mutsizp = dummy$mutsizp
        
        dum = gibbs_admkr_erro(xh = x, inicost = inicost, k = k, errorsizp = errorsizp, errorprob = errorprob, data_x = data_x, data_y = data_y)
        x = dum$x
        inicost = dum$cost
        errorsizpwarm[k] = errorsizp = dum$errorsizp
    }
    return(list(x = x, cost = inicost, mutsizp = mutsizpwarm[warm], errorsizp = errorsizpwarm[warm]))
}

