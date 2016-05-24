
aapc = function(fit, type="HAZ_AC(CS)", interval=5) {

	stopifnot(type %in% c("HAZ_AC(CS)", "HAZ_APC(CS)", "HAZ_APC(HR)"));
	getSeAapc = function(fit, fnAapc) {
		beta = fit$coefficient[, 1];
		delta = fit$coefficient[, 2] * 1e-4;
		est = fnAapc(fit, beta);
		jacob = matrix(NA, length(beta), length(est));
		for (i in 1:length(beta)) {
			beta1 = beta;
			beta1[i] = beta[i] + delta[i];
			est1 = fnAapc(fit, beta1);
			jacob[i, ] = (est1 - est) / delta[i];
		}
		varAapc = t(jacob) %*% fit$covariance %*% jacob;
		if (length(est) == 1) seAapc = sqrt(varAapc)
		else seAapc = sqrt(diag(varAapc));
		result = fit$apc;
		result$estimate = est;
		result$std.error = seAapc;
		result = data.frame(result);
		rownames(result) = NULL;
		return(result);
	}

	fup = interval;
	getAapc = function(fit, beta) {
		epsilon = 1e-5;
		getAapc1 = function(years) {
			pred1 = fit$Predict(years + epsilon, fup, beta_input=beta);
			pred2 = fit$Predict(years - epsilon, fup, beta_input=beta);
			#deriv = (log(pred1$pred_cum) - log(pred2$pred_cum)) / epsilon / 2;
			deriv = ((pred1$pred_cum) - (pred2$pred_cum)) / epsilon / 2;
			return(deriv);
		}
		nSeg = dim(fit$apc)[1];
		newApc = rep(NA, nSeg);
		for (i in 1:nSeg) {
			t0 = fit$apc$start.year[i];
			t1 = fit$apc$end.year[i];
			newApc[i] = mean(getAapc1(t0:t1));
		}
		return(newApc);
	}

	getAapc_log = function(fit, beta) {
		epsilon = 1e-5;
		getAapc1 = function(years) {
			pred1 = fit$Predict(years + epsilon, fup, beta_input=beta);
			pred2 = fit$Predict(years - epsilon, fup, beta_input=beta);
			deriv = (log(pred1$pred_cum) - log(pred2$pred_cum)) / epsilon / 2;
			#deriv = ((pred1$pred_cum) - (pred2$pred_cum)) / epsilon / 2;
			return(deriv);
		}
		nSeg = dim(fit$apc)[1];
		newApc = rep(NA, nSeg);
		for (i in 1:nSeg) {
			t0 = fit$apc$start.year[i];
			t1 = fit$apc$end.year[i];
			newApc[i] = mean(getAapc1(t0:t1));
		}
		return(newApc);
	}

	getAapc_hr = function(fit, beta) {
		nJP = length(fit$jp);
		result = rep(NA, nJP + 1);
		result[1] = beta[nJP + 1];
		if (nJP > 0) for (i in 1:nJP) result[i + 1] = result[i] + beta[i]
		return(exp(result) - 1);
	}

	result = NULL;
	if (type == "HAZ_AC(CS)") result = getSeAapc(fit, getAapc)
	else if (type == "HAZ_APC(CS)") result = getSeAapc(fit, getAapc_log)
	else if (type == "HAZ_APC(HR)") result = getSeAapc(fit, getAapc_hr)
	return(result);
}
