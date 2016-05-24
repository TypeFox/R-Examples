BACprior.lm = function(Y, X, U, omega = c(1, 1.1, 1.3, 1.6, 2, 5, 10, 30, 50, 100, Inf), maxmodels = 150, cutoff = 0.0001, return.best = FALSE)
{
	na.fail(cbind(Y, X, U));
	n = length(X);
	ncov = ncol(U);
	resultsX = summary(regsubsets(y = X, x = U, nbest = maxmodels, really.big = T, nvmax = ncov));
	MLx = exp(-resultsX$bic/2 + min(resultsX$bic)/2);
	MLx = MLx/sum(MLx);
	modelsX = resultsX$which[,-1]; #The null model is not considered

	resultsYa = regsubsets(y = Y, x = cbind(X,U), force.in = 1, nbest = maxmodels, really.big = T, nvmax = ncov + 1);
	resultsY = summary(resultsYa);
	MLy = exp(-resultsY$bic/2 + min(resultsY$bic)/2); 
	MLy = MLy/sum(MLy);
	modelsY = resultsY$which[,-c(1,2)]; #The null model is not considered
	inX_notinY = matrix(apply(modelsY, 1, function(x){apply(modelsX, 1, function(y){sum(y > x)})}), ncol = nrow(modelsX), byrow = T);
		#each column indicate the number of covariates in the different X models that are not in the given Y model
	Nothers = ncov - inX_notinY;
	l = length(omega);
	MF = matrix(0, nrow = nrow(modelsY), ncol = l);
	for(i in 1:l)
	{
		if(omega[i] == Inf)
		{
			a = 1/3;
			b = 0;
		}
		else
		{
			a = omega[i]/(3*omega[i]+1);
			b = 1/(3*omega[i]+1);
		}

		MF[,i] = (a**Nothers*b**inX_notinY)%*%MLx;
	}
	unnormprob = apply(MF, 2, "*", MLy);
	normprob = apply(unnormprob, 2, function(x){x/sum(x)});
	normprob = normprob*(normprob>cutoff);
	usedmodels = which(rowSums(normprob) > 0);
	usedmodels;
	est = coef(resultsYa, id = usedmodels);
	VarEst = vcov(resultsYa, id = usedmodels);
	betas = numeric(nrow(modelsY));
	variances = numeric(nrow(modelsY));
	if(length(usedmodels) == 1)
	{
		betas[usedmodels] = est[2];
		variances[usedmodels] = VarEst[2,2];
	}
	else
	{
		k = 1;
		for(i in usedmodels)
		{
			betas[i] = est[[k]][2];
			variances[i] = VarEst[[k]][2,2];
			k = k + 1;
		}
	}
	beta = apply(normprob, 2, function(x){sum(x*betas)});
	variance1 = apply(normprob, 2, function(x){sum(x*(variances + betas**2))});
	variance2 = variance1 - beta**2;
	results = cbind(omega, beta, sqrt(variance2));
	rownames(results) = NULL;
	colnames(results) = c("omega", "Posterior mean", "Standard deviation");
	if(return.best == FALSE)
	{
		return(list(results = results));
	}
	else
	{
		colnames(modelsY) = colnames(U);
		colnames(normprob) = omega;
		return(list(results = results, best.models = modelsY[usedmodels,], posterior.prob = normprob[usedmodels,]));
	}	
}