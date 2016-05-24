warmup_gaussian = function (x, inicost, mutsizp, warm = 100, prob = 0.234, data_x, data_y, prior_p = 2, prior_st = 1)
{
    mutsizpwarm = vector(, warm)
    for (k in 1:warm)
    {
        dummy = np_gibbs(xh = x, inicost = inicost, k = k, mutsizp = mutsizp, 
            prob = prob, data_x = data_x, data_y = data_y, prior_p = prior_p,
			prior_st = prior_st)
        x = dummy$x
        inicost = dummy$cost
        sigma2 = dummy$sigma2
        mutsizpwarm[k] = mutsizp = dummy$mutsizp
    }
    return(list(x = x, sigma2 = sigma2, cost = inicost, mutsizplast = mutsizpwarm[warm], mutsizp = mutsizpwarm))
}
