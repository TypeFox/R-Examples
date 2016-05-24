mcmcrecord_gaussian = function (x, inicost, mutsizp, warm = 100, M = 100, prob = 0.234, 
    num_batch = 10, step = 10, data_x, data_y, xm, alpha = 0.05, prior_p = 2, 
    prior_st = 1, mlike = c("Chib", "Geweke", "LaplaceMetropolis", "all")) 
{
	mlike = match.arg(mlike)
    dm = ncol(data_x)
    data_num = nrow(data_x)
    hn = data_num^(-1/(4 + dm))
    size_batch = M/num_batch
    sum_h = total_sd = sif = vector(, (dm + 2))
    ave_h2 = h2 = vector(, (dm + 1))
    batch_h = matrix(, num_batch, (dm + 2))
    for (i in 1:(dm + 2))
    {
        sum_h[i] = 0
        for(j in 1:num_batch)
        {
            batch_h[j, i] = 0
        }
        total_sd[i] = 0
    }
    for (i in 1:(dm + 1))
    {
        ave_h2[i] = 0
    }
    sigma2record = costrecord = acceptance = vector(, M)
    xrecord = matrix(, M, dm)
    cpost = matrix(, M/step, dm + 1)
	gpost = matrix(, M, dm + 1)
    for (k in 1:M) {
        dummy = np_gibbs(x, inicost, k + warm, mutsizp, prob = prob, 
            data_x = data_x, data_y = data_y, prior_p = prior_p,
			prior_st = prior_st)
        xrecord[k, ] = x = dummy$x
        sigma2record[k] = dummy$sigma2
        costrecord[k] = inicost = dummy$cost
        acceptance[k] = dummy$accept_h
        sn = ceiling(k/size_batch)
        index = ceiling(k/step)
        sum_h[dm + 2] = sum_h[dm + 2] + costrecord[k]
        batch_h[sn, dm + 2] = batch_h[sn, dm + 2] + costrecord[k]
        total_sd[dm + 2] = total_sd[dm + 2] + (costrecord[k])^2
        sum_h[dm + 1] = sum_h[dm + 1] + sqrt(sigma2record[k])
        ave_h2[dm + 1] = ave_h2[dm + 1] + sigma2record[k]
        batch_h[sn, dm + 1] = batch_h[sn, dm + 1] + sqrt(sigma2record[k])
        total_sd[dm + 1] = total_sd[dm + 1] + sigma2record[k]
        cpost[index, dm + 1] = gpost[k, 1] = sigma2record[k]	 		
        for (j in 1:dm) {
            temp = exp(xrecord[k, j])
            sum_h[j] = sum_h[j] + sqrt(temp)
            ave_h2[j] = ave_h2[j] + temp
            batch_h[sn, j] = batch_h[sn, j] + sqrt(temp)
            total_sd[j] = total_sd[j] + temp
            cpost[index, j] = temp
			gpost[k, j+1] = temp
        }
    }
    std_h = vector(, (dm + 2))
    for (i in 1:(dm + 2))
    {
        std_h[i] = 0
        sum_h[i] = sum_h[i]/M
        for (j in 1:num_batch)
        {
            std_h[i] = std_h[i] + (batch_h[j, i]/size_batch - sum_h[i])^2
        }
        std_h[i] = sqrt(std_h[i]/(num_batch * num_batch - num_batch))
        total_sd[i] = sqrt(M * (total_sd[i]/M - sum_h[i] * sum_h[i])/(M - 1))
        sif[i] = (std_h[i]/total_sd[i])^2 * M
    }
    h2 = ave_h2/M
	if(mlike == "Chib")
	{
		logmargin = loglikelihood_gaussian(h2, data_x, data_y)
		mlikeres = logmargin + logpriors_gaussian(h2, data_x, prior_p = prior_p, prior_st = prior_st) - logdensity_gaussian(h2, cpost)
	}
	if(mlike == "Geweke")
	{
		mlikeres = cov_chol(gpost, data_x, data_y, alpha = alpha, prior_p = prior_p, prior_st = prior_st)
	}
	if(mlike == "LaplaceMetropolis")
	{
		mlikeres = LaplaceMetropolis_gaussian(gpost, data = data_x, data_y = data_y, prior_p = prior_p, prior_st = prior_st)
	}
	else
	{
		logmargin = loglikelihood_gaussian(h2, data_x, data_y)
		Chib = logmargin + logpriors_gaussian(h2, data_x, prior_p = prior_p, prior_st = prior_st) - logdensity_gaussian(h2, cpost)
		Geweke = cov_chol(gpost, data_x, data_y, alpha = alpha, prior_p = prior_p, prior_st = prior_st)
		Laplace = LaplaceMetropolis_gaussian(gpost, data = data_x, data_y = data_y, prior_p = prior_p, prior_st = prior_st)
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
    colnames(sum_h)[dm+1] = colnames(sif)[dm+1] = "sigma"
    colnames(sum_h)[dm+2] = colnames(sif)[dm+2] = "costvalue"
    colnames(h2)[dm+1] = colnames(cpost)[dm+1] = colnames(gpost)[dm+1] = "sigma^2"
    if (missing(xm))
    {
        return(list(sum_h = sum_h, h2 = h2, sif = sif, cpost = cpost, gpost = gpost, 
				accept = mean(acceptance), marginalike = mlikeres))
    }
    else
    {
        ker = kern(sum_h, data_x = data_x, data_y = data_y, xm = xm)
        R2 = ker$r2
        MSE = ker$mse
        return(list(sum_h = sum_h, h2 = h2, sif = sif, cpost = cpost, gpost = gpost, 
				accept = mean(acceptance), marginalike = mlikeres, R2 = R2, MSE = MSE))
    }
}
