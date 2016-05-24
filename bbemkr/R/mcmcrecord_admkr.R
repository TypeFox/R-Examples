mcmcrecord_admkr = function(x, inicost, mutsizp, errorsizp, warm = 100, M = 100, prob = 0.234, errorprob = 0.44, num_batch = 10, step = 10, data_x, data_y,
                       xm, alpha = 0.05, mlike = c("Chib", "Geweke", "LaplaceMetropolis", "all"))
{
    mlike = match.arg(mlike)
    dm = ncol(data_x)
    data_num = nrow(data_x)
    hn = data_num^(-1/(4 + dm))
    bn = data_num^(-3/(2 * dm + 11))
    size_batch = M/num_batch
    len = M/step
    sum_h = total_sd = sif = vector(, dm + 2)
    ave_h2 = h2 = vector(, (dm + 1))
    batch_h = matrix(, num_batch, (dm + 2))
    for(i in 1:(dm+2))
    {
        sum_h[i] = 0
        for(j in 1:num_batch)
        {
            batch_h[j,i] = 0
        }
        total_sd[i] = 0
    }
    for(i in 1:(dm+1))
    {
        ave_h2[i] = 0
    }
    costrecord = acceptance_nw = acceptance_erro = vector(,M)
    xprecord = matrix(, M, (dm + 1))
    cpost = matrix(, M/step, dm + 1)
    gpost = matrix(, M, dm+1)
    for(k in 1:M)
    {
        dummy = gibbs_admkr_nw(xh = x, inicost = inicost, k = k, mutsizp = mutsizp, prob = prob, data_x = data_x, data_y = data_y)
        x = dummy$x
        inicost = dummy$cost
        acceptance_nw[k] = dummy$accept_h
        mutsizp = dummy$mutsizp

        dum = gibbs_admkr_erro(xh = x, inicost = inicost, k = k, errorsizp = errorsizp, errorprob = errorprob, data_x = data_x, data_y = data_y)
        xprecord[k,] = x = dum$x
        costrecord[k] = inicost = dum$cost
        errorsizp = dum$errorsizp
        acceptance_erro[k] = dum$accept_erro

        sn = ceiling(k/size_batch)
        index = ceiling(k/step)
        
        sum_h[dm + 2] = sum_h[dm + 2] + costrecord[k]
        batch_h[sn, dm + 2] = batch_h[sn, dm + 2] + costrecord[k]
        total_sd[dm + 2] = total_sd[dm + 2] + (costrecord[k])^2
        for(j in 1:(dm + 1))
        {
            temp = exp(xprecord[k,j])
            sum_h[j] = sum_h[j] + sqrt(temp)
            ave_h2[j] = ave_h2[j] + temp
            batch_h[sn,j] = batch_h[sn,j] + sqrt(temp)
            total_sd[j] = total_sd[j] + temp
            cpost[index,j] = temp
            gpost[k,j] = temp
        }
    }
    std_h = vector(, dm + 2)
    for(i in 1:(dm+2))
    {
        std_h[i] = 0
        sum_h[i] = sum_h[i]/M
        for(j in 1:num_batch)
        {
            std_h[i] = std_h[i] + (batch_h[j,i]/size_batch - sum_h[i])^2
        }
        std_h[i] = sqrt(std_h[i]/(num_batch^2 - num_batch))
        total_sd[i] = sqrt(M * (total_sd[i]/M - sum_h[i]^2)/(M - 1))
        sif[i] = (std_h[i]/total_sd[i])^2 * M
    }
    h2 = ave_h2/M
    if (mlike == "Chib") {
        logmargin = loglikelihood_admkr(h2, data_x = data_x, data_y = data_y)
        mlikeres = logmargin + logpriors_admkr(h2, data_x = data_x) - logdensity_admkr(h2, cpost)
    }
    if (mlike == "Geweke") {
        mlikeres = cov_chol_admkr(gpost, alpha = alpha)
    }
    if (mlike == "LaplaceMetropolis") {
        mlikeres = LaplaceMetropolis_admkr(gpost, data_x = data_x)
    }
    else {
        logmargin = loglikelihood_admkr(h2, data_x = data_x, data_y = data_y)
        Chib = logmargin + logpriors_admkr(h2, data_x = data_x) - logdensity_admkr(h2, cpost)
        Geweke = cov_chol_admkr(gpost, alpha = alpha, data_x = data_x, data_y = data_y)
        Laplace = LaplaceMetropolis_admkr(gpost, data_x = data_x, data_y = data_y)
        mlikeres = matrix(c(Chib, Geweke, Laplace), ncol = 3)
        colnames(mlikeres) = c("Chib", "Geweke", "Laplace")
        rownames(mlikeres) = ""
    }
    # assigning colnames to variables
    sum_h = data.frame(matrix(sum_h, nrow=1))
    h2 = data.frame(matrix(h2, nrow=1))
    sif = data.frame(matrix(sif, nrow=1))
    cpost = data.frame(cpost)
    gpost = data.frame(gpost)
    for(i in 1:dm)
    {
        colnames(sum_h)[i] = colnames(sif)[i] = paste("h",i,sep="_")
        colnames(h2)[i] = colnames(cpost)[i] = colnames(gpost)[i] = paste(paste("h",i,sep="_"),2,sep="^")
    }
    colnames(sum_h)[dm+1] = colnames(sif)[dm+1] = "b"
    colnames(sum_h)[dm+2] = colnames(sif)[dm+2] = "costvalue"
    colnames(h2)[dm+1] = colnames(cpost)[dm+1] = colnames(gpost)[dm+1] = "b^2"
    if (missing(xm)) {
        return(list(sum_h = sum_h, h2 = h2, sif = sif, cpost = cpost,
            gpost = gpost, accept_nw = mean(acceptance_nw), accept_erro = mean(acceptance_erro), marginalike = mlikeres))
    }
    else {
        ker = kern(sum_h, data_x = data_x, data_y = data_y, xm = xm)
        R2 = ker$r2
        MSE = ker$mse
        return(list(sum_h = sum_h, h2 = h2, sif = sif, cpost = cpost,
        gpost = gpost, accept_nw = mean(acceptance_nw), accept_erro = mean(acceptance_erro), marginalike = mlikeres,
            R2 = R2, MSE = MSE))
    }
}
