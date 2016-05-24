MFDM <- function(mort_female, mort_male, mort_ave, percent_1 = 0.95, percent_2 = 0.95, fh, level = 80, alpha = 0.2, 
				MCMCiter = 100, fmethod = c("auto_arima", "ets"), BC = c(FALSE, TRUE), lambda) 
{
   	fmethod = match.arg(fmethod)
	if(BC == TRUE)
	{
		mort_female = BoxCox(mort_female, lambda)
		mort_male = BoxCox(mort_male, lambda)
		mort_ave = BoxCox(mort_ave, lambda)
	}
	
	# mean
	
	mort_femalemean = rowMeans(mort_female)
	mort_malemean   = rowMeans(mort_male)
	mort_avemean    = rowMeans(mort_ave)
	
	# de-mean
	
	mort_avedemean = mort_ave - matrix(rep(mort_avemean, ncol(mort_female)), ncol = ncol(mort_female))
	mort_femaledemean = mort_female - matrix(rep(mort_femalemean, ncol(mort_female)), ncol = ncol(mort_female))
	mort_maledemean = mort_male - matrix(rep(mort_malemean, ncol(mort_female)), ncol = ncol(mort_female))
	
	# svd
	
	W_1 = scale(t(mort_female), scale = FALSE)
	W_2 = scale(t(mort_male), scale = FALSE)
	mort_avesvd = svd(t(mort_avedemean))
	ncomp_1 = which.min(abs(cumsum(mort_avesvd$d^2/sum(mort_avesvd$d^2)) - percent_1))
	first_percent = (mort_avesvd$d[1])^2/sum(mort_avesvd$d^2)
	basis_ave = mort_avesvd$v[,1:ncomp_1]
	score_ave = t(mort_avedemean) %*% basis_ave
	decomp_ave = basis_ave %*% t(score_ave)
	
	mort_femalediff = t(W_1) - decomp_ave
	mort_malediff = t(W_2) - decomp_ave
	
	# svd (2nd time) (female)
	
	mort_femalesvd = svd(t(mort_femalediff))
	ncomp_female = which.min(abs(cumsum(mort_femalesvd$d^2/sum(mort_femalesvd$d^2)) - percent_2))
	female_percent = (mort_femalesvd$d[1])^2/sum(mort_femalesvd$d^2)
	basis_female = mort_femalesvd$v[,1:ncomp_female]

	mort_malesvd = svd(t(mort_malediff))
	ncomp_male = which.min(abs(cumsum(mort_malesvd$d^2/sum(mort_malesvd$d^2)) - percent_2))
    male_percent = (mort_malesvd$d[1])^2/sum(mort_malesvd$d^2)
	basis_male = mort_malesvd$v[,1:ncomp_male]

	psi_1 = as.matrix(basis_ave)
	psi_2 = as.matrix(basis_female)
	psi_3 = as.matrix(basis_male)

	N_subj = ncol(mort_female)
	N_obs  = nrow(mort_female)
	dim_space_b = ncomp_1
	dim_space_w = ncomp_female
	dim_space_f = ncomp_male
	data_A = list("N_subj", "N_obs", "dim_space_b", "dim_space_w", "dim_space_f", "W_1", "W_2", "psi_1", "psi_2", "psi_3")

	inits_A = function()
	{
		list(ll_w = c(rep(1, ncomp_female)), ll_f = c(rep(1, ncomp_male)), ll_b = c(rep(1, ncomp_1)), 
			taueps_1 = 1, taueps_2 = 1)
	}
		
	zi = matrix(, N_subj, dim_space_w)
	xi = matrix(, N_subj, dim_space_b)
	fi = matrix(, N_subj, dim_space_f)
	ll_b = vector(, dim_space_b)
	ll_w = vector(, dim_space_w)
	ll_f = vector(, dim_space_f)	
	
	inprod = function(a, b)
	{
		return(matrix(a, nrow = 1) %*% matrix(b, ncol = 1))
	}
		
	popmodel = function()
	{
		for(i in 1:N_subj)
		{
			for(t in 1:N_obs)
			{
				W_1[i,t] ~ dnorm(m_1[i,t], taueps_1)
				W_2[i,t] ~ dnorm(m_2[i,t], taueps_2)

				m_1[i,t] <- X[i,t] + U_1[i,t]
				m_2[i,t] <- X[i,t] + U_2[i,t]

				X[i,t] <- inprod(xi[i,], psi_1[t,])
				U_1[i,t] <- inprod(zi[i,], psi_2[t,])
				U_2[i,t] <- inprod(fi[i,], psi_3[t,])				
			}
			for(k in 1:dim_space_b)
			{
				xi[i,k] ~ dnorm(0.0, ll_b[k])
			}
			for(l in 1:dim_space_w)
			{
				zi[i,l]  ~ dnorm(0.0, ll_w[l])
			}
			for(j in 1:dim_space_f)
			{	
				fi[i,j] ~ dnorm(0.0, ll_f[j])
			}
		}	
		for(k in 1:dim_space_b)
		{
			ll_b[k] ~ dgamma(1.0E-3, 1.0E-3)
			lambda_b[k] <- 1/ll_b[k]
		}
		for(l in 1:dim_space_w)
		{
			ll_w[l] ~ dgamma(1.0E-3, 1.0E-3)
			lambda_w[l] <- 1/ll_w[l]
		}
		for(j in 1:dim_space_f)
		{
			ll_f[j] ~ dgamma(1.0E-3, 1.0E-3)
			lambda_f[j] <- 1/ll_f[j]
		}
		
		# prior

		taueps_1 ~ dgamma(1.0E-3, 1.0E-3)
		taueps_2 ~ dgamma(1.0E-3, 1.0E-3)
	}	
	if(requireNamespace("R2jags", quietly = TRUE)) 
	{
		bugs_chose = R2jags::jags(data = data_A, inits = inits_A, model.file = popmodel,
		                    parameters.to.save = c("xi", "zi", "fi", "taueps_1", "taueps_2", "lambda_b", "lambda_f", "lambda_w"),
		                    n.chains = 1, n.burnin = 5000, n.iter = 6000, DIC = TRUE)
	}
	else
	{
		stop("Please install JAGS")
	}
	mortality_female_coda = bugs_chose$BUGSoutput$sims.list
		
    taueps_female = mortality_female_coda$taueps_1
	taueps_male   = mortality_female_coda$taueps_2
		
    lambda_coda_common = mortality_female_coda$lambda_b
	lambda_coda_male   = mortality_female_coda$lambda_f
	lambda_coda_female = mortality_female_coda$lambda_w

	score_common = array(,c(MCMCiter,ncol(mort_female),ncomp_1))
	score_female = array(,c(MCMCiter,ncol(mort_female),ncomp_female))
	score_male = array(,c(MCMCiter,ncol(mort_male),ncomp_male))
	
	if(ncomp_1 == 1)
    {
       for(i in 1:MCMCiter)
       {
           score_common[i,,1] = mortality_female_coda$xi[i,]
       }
    }
    else
    {
       for(i in 1:MCMCiter)
       {
           score_common[i,,] = mortality_female_coda$xi[i,,]
       }
	}
        
    if(ncomp_female == 1)
    {
       for(i in 1:MCMCiter)
       {
           score_female[i,,1] = mortality_female_coda$zi[i,]
       }
    }
    else
    {
       for(i in 1:MCMCiter)
       {
           score_female[i,,] = mortality_female_coda$zi[i,,]
       }
	}
        
    if(ncomp_male == 1)
    {
       for(i in 1:MCMCiter)
       {
           score_male[i,,1] = mortality_female_coda$fi[i,]
       }
    }
    else
    {
       for(i in 1:MCMCiter)
       {
           score_male[i,,] = mortality_female_coda$fi[i,,]
       }
    }
	
	qconf = qnorm(.5 + level/200)
	score_common_fore = score_common_varfcast = array(, dim = c(MCMCiter, fh, ncomp_1))		
	for(i in 1:MCMCiter)
	{
		if(fmethod == "auto_arima")
		{
			for(j in 1:ncomp_1)
			{
				dum = forecast(auto.arima(score_common[i,,j]), h = fh, level = level)
				score_common_varfcast[i,,j] = ((dum$upper - dum$lower)/(2*qconf))^2
				score_common_fore[i,,j] = dum$mean
			}
		}
		if(fmethod == "ets")
		{
			for(j in 1:ncomp_1)
			{
				dum = forecast(ets(score_common[i,,j]), h = fh, level = level)
				score_common_varfcast[i,,j] = ((dum$upper - dum$lower)/(2*qconf))^2
				score_common_fore[i,,j] = dum$mean
			}			
		}
	}

	score_female_fore = score_female_varfcast = array(,dim=c(MCMCiter, fh, ncomp_female))
	for(i in 1:MCMCiter)
	{
		if(fmethod == "auto_arima")
		{
			for(j in 1:ncomp_female)
			{
				dum = forecast(auto.arima(score_female[i,,j]), h = fh, level = level)
				score_female_varfcast[i,,j] = ((dum$upper - dum$lower)/(2*qconf))^2
				score_female_fore[i,,j] = dum$mean		
			}
		}
		if(fmethod == "ets")
		{
			for(j in 1:ncomp_female)
			{
				dum = forecast(ets(score_female[i,,j]), h = fh, level = level)
				score_female_varfcast[i,,j] = ((dum$upper - dum$lower)/(2*qconf))^2
				score_female_fore[i,,j] = dum$mean		
			}			
		}
	}
		
	score_male_fore = score_male_varfcast = array(,dim = c(MCMCiter, fh, ncomp_male))
	for(i in 1:MCMCiter)
	{
		if(fmethod == "auto_arima")
		{
			for(j in 1:ncomp_male)
			{	
				dum = forecast(auto.arima(score_male[i,,j]), h = fh, level = level)
				score_male_varfcast[i,,j] = ((dum$upper - dum$lower)/(2*qconf))^2
				score_male_fore[i,,j] = dum$mean
			}
		}
		if(fmethod == "ets")
		{
			for(j in 1:ncomp_male)
			{	
				dum = forecast(ets(score_male[i,,j]), h = fh, level = level)
				score_male_varfcast[i,,j] = ((dum$upper - dum$lower)/(2*qconf))^2
				score_male_fore[i,,j] = dum$mean
			}
		}		
	}
	
	forescore_common = array(, dim = c(MCMCiter, fh, ncomp_1))
	for(i in 1:MCMCiter)
	{
		for(j in 1:ncomp_1)
		{	
			for(h in 1:fh)
			{
				forescore_common[i,h,j] = rnorm(1, mean = score_common_fore[i,h,j], 
												sd = sqrt(score_common_varfcast[i,h,j]))
			}
		}			
	}
			
	forescore_female = array(, dim = c(MCMCiter, fh, ncomp_female))
	for(i in 1:MCMCiter)
	{
		for(j in 1:ncomp_female)
		{
			for(h in 1:fh)
			{
				forescore_female[i,h,j] = rnorm(1, mean = score_female_fore[i,h,j], 
												sd = sqrt(score_female_varfcast[i,h,j]))
			}
		}
	}

	forescore_male = array(, dim = c(MCMCiter, fh, ncomp_male))
	for(i in 1:MCMCiter)
	{
		for(j in 1:ncomp_male)
		{
			for(h in 1:fh)
			{
				forescore_male[i,h,j] = rnorm(1, mean = score_male_fore[i,h,j],
												sd = sqrt(score_male_varfcast[i,h,j]))
			}
		}
	}
	
	ave_boot_fore = array(, dim = c(MCMCiter, fh, nrow(mort_female)))
	for(i in 1:MCMCiter)
	{
		for(h in 1:fh)
		{
			ave_boot_fore[i,h,] = t((forescore_common[i,h,] %*% t(basis_ave)))
		}
	}
			
	female_boot_fore = array(, dim = c(MCMCiter, fh, nrow(mort_female)))
	for(i in 1:MCMCiter)
	{
		for(h in 1:fh)
		{
			female_boot_fore[i,h,] = t((forescore_female[i,h,] %*% t(basis_female))) 
		}
	}
			
	male_boot_fore = array(, dim = c(MCMCiter, fh, nrow(mort_female)))
	for(i in 1:MCMCiter)
	{
		for(h in 1:fh)
		{
			male_boot_fore[i,h,] = t((forescore_male[i,h,] %*% t(basis_male)))
		}
	}	
	
	female_fore = array(, dim = c(MCMCiter, fh, nrow(mort_female)))
	for(i in 1:MCMCiter)
	{
		for(h in 1:fh)
		{
			female_fore[i,h,] = mort_femalemean + ave_boot_fore[i,h,] + female_boot_fore[i,h,]
		}
	}

	male_fore = array(, dim = c(MCMCiter, fh, nrow(mort_male)))
	for(i in 1:MCMCiter)
	{
		for(h in 1:fh)
		{
			male_fore[i,h,] = mort_malemean + ave_boot_fore[i,h,] + male_boot_fore[i,h,]
		}
	}	
	
	female_fore_variance = array(, dim = c(MCMCiter, fh, nrow(mort_female)))
	for(i in 1:MCMCiter)
	{
		for(h in 1:fh)
		{
			for(j in 1:nrow(mort_female))
			{
				female_fore_variance[i,h,j] = rnorm(1, mean = female_fore[i,h,j], 
													sd = sqrt(1/taueps_female[i]))			
			}
		}	
	}
			
	male_fore_variance = array(, dim = c(MCMCiter, fh, nrow(mort_male)))
	for(i in 1:MCMCiter)
	{
		for(h in 1:fh)
		{
			for(j in 1:nrow(mort_male))
			{
				male_fore_variance[i,h,j] = rnorm(1, mean = male_fore[i,h,j],
													sd = sqrt(1/taueps_male[i]))		
			}
		}
	}
	
	if(BC == TRUE)
	{
		mort_female_fore = InvBoxCox(female_fore_variance, lambda)
		mort_male_fore   = InvBoxCox(male_fore_variance, lambda)
	}
	else
	{
		mort_female_fore = female_fore_variance
		mort_male_fore   = male_fore_variance
	}
	return(list(first_percent = first_percent, female_percent = female_percent, male_percent = male_percent,
                    mort_female_fore = mort_female_fore, mort_male_fore = mort_male_fore))
}

