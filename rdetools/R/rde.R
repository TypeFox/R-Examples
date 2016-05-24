`rde` <-
function(K, y, est_y = FALSE, alldim = FALSE, est_noise = FALSE, regression = FALSE, nmse = TRUE, dim_rest = 0.5, tcm = TRUE)
{
	if(!tcm)
	{
		# rde by leave-one-out cross-validation (loo-cv)
		return(rde_loocv(K = K, y = y, est_y = est_y, alldim = alldim, est_noise = est_noise, regression = regression, nmse = nmse, dim_rest = dim_rest))
	}
	else
	{
		# rde by fitting a two-component model (tcm)
		return(rde_tcm(K = K, y = y, est_y = est_y, alldim = alldim, est_noise = est_noise, regression = regression, nmse = nmse, dim_rest = dim_rest))
	}
}

