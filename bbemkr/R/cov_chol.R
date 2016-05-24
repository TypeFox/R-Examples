cov_chol <- function(xpost, data_x, data_y, alpha, prior_p, prior_st)
{
    dm = ncol(xpost)
    M = nrow(xpost)
    ave = temv = vector(,dm)
    b = cov = matrix(,dm,dm)
    fhat = tail = vector(,M)

    ave = colMeans(xpost)
    for(i in 1:dm)
    {
        for(j in 1:i)
        {
            tem = 0
            for(k in 1:M)
            {
                tem = tem + (xpost[k,i] - ave[i]) * (xpost[k,j] - ave[j])
            }
            cov[i,j] = tem/M
        }
    }
    for(i in 1:(dm-1))
    {
        for(j in (i+1):dm)
        {
            cov[i,j] = cov[j,i]
        }
    }
    b = solve(chol(cov))
    det = det(b)
    fmean = iave = 0
    for(i in 1:M)
    {
        temv[1] = b[1,1] * (xpost[i,1] - ave[1])
        for(j in 2:dm)
        {
            tem = 0
            for(k in 1:j)
            {
                tem = tem + b[j,k] * (xpost[i,k] - ave[k])
            }
            temv[j] = tem
        }
        tem = 0
        for(j in 1:dm)
        {
            tem = tem + temv[j]^2
        }
        if(tem < qchisq(1-alpha, dm))
        {
            fhat[i] = -0.5 * dm * log(2.0 * pi) + log(det) - 0.5 * tem - log(1.0 - alpha) -
                        margin_prior_gaussian(xpost[i,], data_x = data_x, prior_p = prior_p, prior_st = prior_st) -
                        margin_like_gaussian(xpost[i,], data_x = data_x, data_y = data_y)
            tail[i] = 1
        }
        else
        {
            fhat[i] = 0
            iave = iave + 1
            tail[i] = 0
        }
        fmean = fmean + fhat[i]/M
    }
    fmean = fmean/(M - iave) * M
    tem = 0
    for(i in 1:M)
    {
        tem = tem + exp(tail[i] * (fhat[i] - fmean)) * tail[i]
    }
    tem = tem/M
    tem = log(1.0/tem) - fmean
    return(tem)
}

