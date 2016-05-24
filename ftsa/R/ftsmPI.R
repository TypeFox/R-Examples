ftsmPI <- function (object, B, level, h, fmethod = c("ets", "arima")) 
{
    data = object$y$y
    p = dim(data)[1]
    n = dim(data)[2]
    ncomp = dim(object$basis)[2] - 1
    mdata = apply(data, 1, mean)
    mdata2 = array(rep(as.matrix(mdata), B * h), dim = c(p, B, h))
    sdata = scale(t(data), scale = FALSE)
    load = svd(sdata)$v[, 1:ncomp]
    sco = sdata %*% load
    olivia = matrix(, ncomp, h)
    if(fmethod == "ets")
    {
        for(i in 1:ncomp) 
        {
            olivia[i, ] = forecast(ets(sco[, i]), h = h)$mean
        }
    }
    if(fmethod == "arima")
    {
        for(i in 1:ncomp) 
        {
            olivia[i, ] = forecast(auto.arima(sco[, i]), h = h)$mean
        }
    }
    # in-sample forecast error
    forerr = matrix(, (n - ncomp - h + 1), ncomp)
    for(i in h:(n - ncomp)) 
    {
        k = i + (ncomp - h)
        fore = matrix(, 1, ncomp)
        if(fmethod == "ets") 
        {
            for(j in 1:ncomp) 
            {
                fore[, j] = forecast(ets(sco[1:k, j]), h = h)$mean[h]
            }
        }
        if(fmethod == "arima")
        {
            for(j in 1:ncomp) 
            {
                fore[, j] = forecast(auto.arima(sco[1:k, j]), h = h)$mean[h]
            }
        }
        forerr[i-h+1, ] = sco[k + h, ] - fore
    }
    # bootstrap residuals (sampling with replacement)
    resi = t(sdata) - load %*% t(sco)
    q = array(, dim = c(p, B, h))
    for(j in 1:h) 
    {
        for(i in 1:p) 
        {
            q[i, , j] = sample(resi[i, ], size = B, replace = TRUE)
        }
    }
    # bootstrap forecast error of principal component scores (sampling with replacement)
    ny = array(, dim = c(ncomp, B, h))
    for(j in 1:h) 
    {
        for(i in 1:ncomp) 
        {
            ny[i, , j] = sample(forerr[, i], size = B, replace = TRUE)
        }
    }
    # add forecasts of principal component scores and bootstrap forecast errors
    oli = array(rep(olivia, B * h), dim = c(ncomp, B, h))
    fo = array(, dim = c(ncomp, B, h))
    for(j in 1:h) 
    {
        for(i in 1:B) 
        {
            fo[, i, j] = oli[, i, j] + ny[, i, j]
        }
    }
    pred = array(, dim = c(p, B, h))
    for(j in 1:h) 
    {
        for(i in 1:B) 
        {
            pred[, i, j] = load %*% fo[, i, j] + mdata2[, i, j] + q[, i, j]
        }
    }
    k1 = k2 = matrix(, p, h)
    for(j in 1:h) 
    {
        for(i in 1:p) 
        {
            k1[i, j] = quantile(pred[i, , j], (100 - level)/200, na.rm = TRUE)
            k2[i, j] = quantile(pred[i, , j], 1 - (100 - level)/200, na.rm = TRUE)
        }
    }
    return(list(bootsamp = pred, lb = k1, ub = k2))
}
