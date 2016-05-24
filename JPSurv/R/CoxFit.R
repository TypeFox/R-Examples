#####################################################################
# The purpose of this program is to estimate and to forecast:
# 1. Proportional hazard relative survival models 
# 2. Survival join point models
# The estimates are computed using the iteratively reweighted least square algorithm.
# The join points are selected using the Bayesian information criterion.
#####################################################################

# CoxFit() is the core function to fit the proportional hazard relative survival models.
# Inputs:
#     X:     the design matrix
#     nAlive: number of patients at risk
#     nDied:  number of people died
#     nLost:  number of peple lost to follow-up
#     expSurv:expected survival

CoxFit = function(X, nAlive, nDied, nLost, expSurv) {
	stopifnot(sum(expSurv <= 0) == 0);
	N = nAlive - nLost / 2;
	Y = N - nDied;
	p = dim(X)[2];

	min_value = 1e-6;
	min_value2 = 1e-15;
	
	# RelativeHazardLinkFun is the link function for the relative hazard model.
	RelativeHazardLinkFun = list();
	RelativeHazardLinkFun$LogLik = function(beta1) {
		xbeta1 = X %*% beta1;
		prob = exp(log(expSurv) * exp(xbeta1));
		nlogl = sum(Y * (log(prob + min_value2) - log((Y+min_value2)/N)) 
				+ (N - Y) * (log(1 - prob + min_value2) - log((N-Y+min_value2)/N)));
		return(nlogl);
	}
	
	RelativeHazardLinkFun$ZW = function(beta) {
		xbeta = X %*% beta;
		prob = exp(log(expSurv) * exp(xbeta));
		prob[prob > 1 - min_value] = 1-min_value;
		prob[prob < min_value] = min_value;
		grad_ = 1 / prob / log(prob);
		Z = (Y / N - prob) * grad_;
		W = prob * (1 - prob) / N * grad_ * grad_;
		
		result = list();
		result$Z = Z;
		result$W = W;
		return(result);
	}
	
	# RelativeSurvivalLinkFun is the link function for the relative survival model.
	RelativeSurvivalLinkFun = list();
	RelativeSurvivalLinkFun$LogLik = function(beta1) {
		xbeta1 = X %*% beta1;
		prob = exp(-exp(xbeta1)) * expSurv;
		nlogl = sum(Y * (log(prob + min_value2) - log((Y+min_value2)/N)) 
				+ (N - Y) * (log(1 - prob + min_value2) - log((N-Y+min_value2)/N)));
		return(nlogl);
	}
	
	RelativeSurvivalLinkFun$ZW = function(beta) {
		xbeta = X %*% beta;
		prob = exp(-exp(xbeta)) * expSurv;
		index = (prob > expSurv - min_value);
		prob[index] = (expSurv - min_value)[index];
		prob[prob < min_value] = min_value;
		grad_ = 1 / prob / log(prob / expSurv);
		Z = (Y / N - prob) * grad_;
		W = prob * (1 - prob) / N * grad_ * grad_;
		
		result = list();
		result$Z = Z;
		result$W = W;
		return(result);
	}

	RelativeSurvivalLinkFun$f = function(xbeta) {
		result = exp(-exp(xbeta));
		return(result);
	}
	
	# GetEstimate() is the function to compute the estimates using the iteratively reweighted least square algorithm.
	GetEstimate = function(linkfun) {
		old_beta = rep(0, p);
		old_ll = linkfun$LogLik(old_beta);
		converged = FALSE;
		lm_fit = NULL;
		for (i in 1:100) {
			ZW = linkfun$ZW(old_beta);
			lm_fit = lm(ZW$Z ~ -1 + X, weights = 1 / ZW$W);
			delta = lm_fit$coefficients;
			scale = 1.0;
			new_beta = NA;
			new_ll = NA;
			for (i in 1:20) {
				new_beta = old_beta + delta * scale;
				new_ll = linkfun$LogLik(new_beta);
				if (new_ll > old_ll) break
				scale = scale / 2;
			}
			if (new_ll < old_ll + 1e-4) {
				if (new_ll > old_ll) {
					converged = TRUE;
				}
				break;
			}
			old_beta = new_beta;
			old_ll = new_ll;
		}
		
		if (i == 100) converged = FALSE;
		
		result = list();
		#result$beta = new_beta;
		result$xbeta = X %*% new_beta;
		#pred = linkfun$f(result$xbeta);
		#result$pred = pred;
		result$ll = new_ll;
		Estimates = new_beta;
		summ_fit = summary(lm_fit);
		Std.Error = summ_fit$coefficients[, "Std. Error"];
		result$coefficients = cbind(Estimates, Std.Error);
		result$covariance = vcov(lm_fit);
		result$aic = 2 * p - 2 * new_ll;
		result$bic = p * log(sum(N)) -2 * new_ll;
		result$converged = converged;
		
		return(result);
	}

	fit = GetEstimate(RelativeSurvivalLinkFun);
	return(fit);
}

# CoxModel_Year() is a wrap-up of the core function CoxFit() when the only covariate is "year of diagnosis".
# Inputs:
#    formu: the forumula in the form ~nAlive+nDied+nLost+expSurv+Interval+year
#    dataSource: the data source
CoxModel_Year = function(formula, data, subset, ...) {
    mf <- match.call(expand.dots=FALSE)
	m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    mt <- attr(mf,"terms")
    dataMatrix = model.frame(mt, mf);
	nAlive = dataMatrix[, 1];
	nDied = dataMatrix[, 2];
	nLost = dataMatrix[, 3];
	ExpSurv = dataMatrix[, 4];
	Interval = dataMatrix[, 5];
	Year = dataMatrix[, 6];

	n = length(Year);
	tempx = log(-log(ExpSurv));
	Interval_ = as.factor(Interval);
	X = model.matrix(~-1+Year+Interval_);
	cox_fit = CoxFit(X, nAlive, nDied, nLost, ExpSurv);
	return(cox_fit);
}

# JoinPointModel_Year() fits the survival joint model when the only covariate is "year of diagnosis".
# The criteria for selecting the join points is Bayesian information criterion.
# Inputs:
#    formu: the forumula in the form ~nAlive+nDied+nLost+expSurv+Interval+year
#    dataSource: the data source
joinpoint = function(formula, data, subset, numJPoints = 0, ...) {
	is.seerstat = FALSE;
	if (length(all.vars(formula)) == 1) {
		varname = all.vars(formula);
		formu = paste("~", "Alive_at_Start + Died + Lost_to_Followup + Expected_Survival_Interval + Interval +", 
			varname, "+ Relative_Survival_Cum", sep = " ");
		formula = as.formula(formu);
		is.seerstat = TRUE;
	}
    mf <- match.call(expand.dots=FALSE)
	m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$formula = formula;
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    mt <- attr(mf,"terms")
    dataMatrix = model.frame(mt, mf);
    remove(subset);
	stopifnot(dim(dataMatrix)[2] == 6 || dim(dataMatrix)[2] == 7);
	nAlive = dataMatrix[, 1];
	nDied = dataMatrix[, 2];
	nLost = dataMatrix[, 3];
	ExpSurv = dataMatrix[, 4];
	Interval = dataMatrix[, 5];
	Year = dataMatrix[, 6];
	RelSurvCum = NULL;
	if (is.seerstat) RelSurvCum = dataMatrix[, "Relative_Survival_Cum"];
	n = length(Year);
	tempx = log(-log(ExpSurv));
	Interval_ = as.factor(Interval);
	X = model.matrix(~-1+Year+Interval_);
	
	# indv_bic() computes the BIC for the given join points.
	indv_bic = function(jpoints) {
		lineSeg = NULL;
		nJP = length(jpoints);
		if (nJP > 0) {
			for (i in 1:nJP) {
				seg = Year - jpoints[i];
				seg[seg < 0] = 0;
				lineSeg = cbind(lineSeg, seg);
			}
			colnames(lineSeg) = paste("jp", 1:nJP, sep = "_");
		}
		
		X1 = cbind(lineSeg, X)
		cox_fit = CoxFit(X1, nAlive, nDied, nLost, ExpSurv);
		
		# predict() computes the predicted cumulative survival rates.
		predict = function(years = NULL, intervals = NULL, beta_input = NULL) {
			if (is.null(years)) years = min(Year):(max(Year) + 10)
			if (is.null(intervals)) intervals = (1:max(Interval))
			if (is.null(beta_input)) beta_input = cox_fit$coefficients[, 1];
			maxInterval = max(intervals);
			stopifnot(maxInterval <= max(Interval));
			
			beta2 = list();
			#cox_fit$beta = cox_fit$coefficients[, 1];
			#intercept = cox_fit$beta[nJP + 1];
			intercept = 0;
			if (nJP > 0) beta2$Seg = beta_input[1:nJP]
			beta2$Year = beta_input[nJP + 1];
			beta2$Interval = c(0, beta_input[(nJP+2):length(beta_input)]);
			beta2$Interval = beta2$Interval + intercept;
			beta2$Interval = beta_input[(nJP+2):length(beta_input)];
			
			surv_int = matrix(NA, length(years), maxInterval);
			surv_cum = surv_int;
			surv_xbeta = surv_int;
			result = data.frame(matrix(0, length(years) * maxInterval, 4));
			names(result) = c("Year", "Interval", "pred_int", "pred_cum");
			for (i in 1:length(years)) {
				for (j in 1:maxInterval) {
					xbeta2 = beta2$Interval[j] + years[i] * beta2$Year;
					if (nJP > 0) for (k in 1:nJP) {
						xseg = years[i] - jpoints[k];
						if (xseg < 0) xseg = 0
						xbeta2 = xbeta2 + beta2$Seg[k] * xseg;
					}
					surv_xbeta[i, j] = xbeta2;
					surv_int[i, j] = exp(-exp(xbeta2));
					if (j == 1) surv_cum[i, j] = surv_int[i, j]
					else surv_cum[i, j] = surv_cum[i, j-1] * surv_int[i, j]
					idx = (i - 1) * maxInterval + j;
					result[idx, 1] = years[i];
					result[idx, 2] = j;
					result[idx, 3] = surv_int[i, j];
					result[idx, 4] = surv_cum[i, j];
				}
			}
			result1 = subset(result, Interval %in% intervals);
			return(result1);
		}
		
		getApc = function() {
			result = matrix(NA, nJP+1, 3);
			result = data.frame(result);
			dimnames(result)[[2]] = c("start.year", "end.year", "estimate");
			#dimnames(result)[[1]] = paste("Seg", 1:(nJP+1), sep = " ");
			rownames(result) = NULL;
			result$estimate[1] = cox_fit$coefficients[nJP+1, 1];
			result$start.year[1] = min(Year);
			result$end.year[nJP + 1] = max(Year);
			if (nJP > 0) {
				for (i in 1:nJP) {
					result$estimate[i + 1] = result$estimate[i] + cox_fit$coefficients[i, 1];
					result$end.year[i] = jpoints[i];
					result$start.year[i + 1] = jpoints[i];
				}
			}
			return(result);
		}
	
		cox_fit$Predict = predict;
		cox_fit$apc = getApc();
		return(cox_fit);
	}
	
	getJP_0 = function() {
		result = indv_bic(NULL);
		result$jp = NULL;
		return(result);
	}
	
	getJP_1 = function() {
		result = list();
		if (max(Year) < min(Year) + 6) return(NA)
		result$jp = NA;
		result$bic = Inf;
		#for (jp1 in (min(Year)+3):(max(Year)-3)) {
		for (jp1 in 5:6) {
			new_bic = indv_bic(jp1);
			if (new_bic$bic < result$bic) {
				result = new_bic;
				result$jp = jp1;
			}
		}
		
		return(result);
	}	
	
	getJP_2 = function() {
		result = list(); 
		if (max(Year) < min(Year) + 9) return(NA)
		result$jp = NA;
		result$bic = Inf;
		#for (jp1 in (min(Year)+3):(max(Year)-6)) 
		#for (jp2 in (jp1 + 3):(max(Year) - 3)) {
		for (jp1 in 13:14) 
		for (jp2 in 17:19) {
			new_bic = indv_bic(c(jp1, jp2));
			if (new_bic$bic < result$bic) {
				result = new_bic;
				result$jp = c(jp1, jp2);
			}
		}
		return(result);
	}	
		
	
	getJP = function(nJP) {
		intervalSize = 3;
		stopifnot(max(Year) > min(Year) + (nJP + 1) * intervalSize);
		jpoints = min(Year) + (1:nJP) * intervalSize;
		endFlag = FALSE;
		nextjp = function(oldjp) {
			oldjp[nJP] = oldjp[nJP] + 1;
			if (nJP >= 2) for (k in nJP:2) {
				if (oldjp[k] > max(Year) - intervalSize * (nJP - k + 1)) {
					oldjp[k - 1] = oldjp[k - 1] + 1;
					oldjp[k:nJP] = oldjp[k - 1] + intervalSize * (1:(nJP - k + 1));
				}
			}
			if (oldjp[1] >= max(Year) - intervalSize * nJP) endFlag <<- TRUE;
			return(oldjp);
		}
		result = list(); 
		result$jp = NA;
		result$bic = Inf;
		nIter = 0;
		jp.debug = FALSE;  # setting the debug option
		while (!endFlag){
			new_bic = indv_bic(jpoints);
			if (new_bic$bic < result$bic) {
				result = new_bic;
				result$jp = jpoints;
			}
			jpoints = nextjp(jpoints);
			nIter = nIter + 1;
			if (jp.debug & nIter > 3) break
		}
		return(result);
	}
	
	getBestJP = function(nJP) {
		bestFit = getJP_0();
		if (nJP > 0) {
			for (k in 1:nJP) {
				cat("Computing estimates when number of join points is ");
				cat(k);
				cat("\n");
				newFit = getJP(k);
				if (newFit$bic < bestFit$bic) bestFit = newFit
			}
		}
		return(bestFit);
	}
	
	cox_fit = getBestJP(numJPoints);

	plot.surv = function() {
		stopifnot(is.seerstat);
		# Get the predicted values
		pred = predict(cox_fit);
		# Plot the predicted values of the 1-, 3-, 5-year survival
		pred1 = subset(pred, Interval == 1);
		plot(pred1$Year, pred1$pred_cum, type = "l", lwd = 3, col = "black", xlab = "Year of diagnosis", ylab = "relative survival", ylim = c(0, 1));
		pred2 = subset(pred1, Year %in% (cox_fit$jp));
		lines(pred2$Year, pred2$pred_cum, type = "p", pch = 20, cex = 1.8, col = "black");
		pred1 = subset(pred, Interval == 3);
		lines(pred1$Year, pred1$pred_cum, type = "l", lwd = 3, col = "darkgreen");
		pred2 = subset(pred1, Year %in% (cox_fit$jp));
		lines(pred2$Year, pred2$pred_cum, type = "p", pch = 20, cex = 1.8, col = "darkgreen");
		pred1 = subset(pred, Interval == 5);
		lines(pred1$Year, pred1$pred_cum, type = "l", lwd = 3, col = "darkorange");
		pred2 = subset(pred1, Year %in% (cox_fit$jp));
		lines(pred2$Year, pred2$pred_cum, type = "p", pch = 20, cex = 1.8, col = "darkorange");
	
		# Plot the observed values of the 1-, 3-, 5-year survival
		lifeDat2 = subset(dataMatrix, Interval == 1);
		lines(lifeDat2[, 6], lifeDat2[, 7], type = "p", col = "black");
		lifeDat2 = subset(dataMatrix, Interval == 3);
		lines(lifeDat2[, 6], lifeDat2[, 7], type = "p", col = "darkgreen");
		lifeDat2 = subset(dataMatrix, Interval == 5);
		lines(lifeDat2[, 6], lifeDat2[, 7], type = "p", col = "darkorange");

		legend("bottomright", c("1-yr est.", "3-yr est.", "5-yr est", "1-yr obs.", "3-yr obs.", "5-yr obs."), 
			lty = c(1, 1, 1, 3, 3, 3), lwd = 3, col = c("black", "darkgreen", "darkorange", "black", "darkgreen", "darkorange"));
	}
	
	pred = cox_fit$Predict();
	cox_fit$predicted = merge(dataMatrix, pred, by.x = c(6, 5), by.y = c("Year", "Interval"));
	cox_fit$plot.surv = plot.surv;
	class(cox_fit) = "joinpoint";
	return(cox_fit);	
}


print.joinpoint = function(x, ...) {
	object = x;
	nJP = length(object$jp);
	cat("Annual percentage changes:\n");
	print(object$apc);
	cat("Coefficient estimates:\n");
	print(object$coefficients[1:(nJP + 1), ]);
	cat("\nInterval estimates:\n");
	print(object$coefficients[-(1:(nJP + 1)), ]);
}

summary.joinpoint = function(object, ...) {
	result = list();
	nJP = length(object$jp);
	result$jp = object$jp;
	result$apc = object$apc;
	result$beta_year = object$coefficients[nJP+1, ];
	result$beta_jp = object$coefficients[1:nJP, 1];
	result$beta_int = object$coefficients[-(1:(nJP + 1)), ];
	return(result);
}

predict.joinpoint = function(object, ...) {
	result = object$Predict(...);
	return(result);
}

plot.joinpoint = function(x, ...) {
	x$plot.surv();
}
