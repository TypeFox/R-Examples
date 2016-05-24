#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

# BDS Test
# Dechert, Scheinkman and LeBaron (1996)
# BDS i.i.d test of 
# Brock and Potter (1993) and de Lima (1996) show that is
# applying the BDS test to the log of the squared standardized residuals
# log(z^2) the bias in the BDS test is almost corrected because the logarithmic
# transformation transforms the GARCH model into a linear additive model which
# satisfies the nuisance free parameter conditions in de Lima (1996) for the BDS
# test (does not apply to more richly parametrized models e.g. with leverage since
# they cannot be cast into linear additive models).

#.bds.test = function(z){
#	fNonlinear::bdsTest(log(z[is.finite(z)]^2), m = 5)
#}

.information.test = function(LLH, nObs, nPars)
{
	AIC  = (-2*LLH)/nObs + 2 * nPars/nObs
	BIC  = (-2*LLH)/nObs + nPars * log(nObs)/nObs
	SIC  = (-2*LLH)/nObs + log((nObs+2*nPars)/nObs)
	HQIC = (-2*LLH)/nObs + (2*nPars*log(log(nObs)))/nObs
	informationTests = list(AIC = AIC, BIC = BIC, SIC = SIC, HQIC = HQIC)
	return(informationTests)
}

# Q-Statistics on Standardized Residuals
.box.test = function(stdresid, p=1, df = 0)
{
	if(any(!is.finite(stdresid))) stdresid[!is.finite(stdresid)]=0
	# p=1 normal case, p=2 squared std. residuals
	# Q-Statistics on Standardized Residuals
	#H0 : No serial correlation ==> Accept H0 when prob. is High [Q < Chisq(lag)]
	box10 = Box.test(stdresid^p, lag = 1, type = "Ljung-Box", fitdf = 0)
	box15 = Box.test(stdresid^p, lag = df+1, type = "Ljung-Box", fitdf = df)
	box20 = Box.test(stdresid^p, lag = df+5, type = "Ljung-Box", fitdf = df)
	LBSR<-matrix(NA,ncol=2,nrow=3)
	LBSR[1:3,1] = c(box10$statistic[[1]],box15$statistic[[1]],box20$statistic[[1]])
	LBSR[1:3,2] = c(box10$p.value[[1]],box15$p.value[[1]],box20$p.value[[1]])
	rownames(LBSR) = c(paste("Lag[1]",sep=""), paste("Lag[p+q+1][",df+1,"]",sep=""), paste("Lag[p+q+5][",df+5,"]",sep=""))
	colnames(LBSR) = c("statistic","p-value")
	return(LBSR)
}

.weightedBoxTest = function(stdresid, p=1, df = 0)
{
	if(any(!is.finite(stdresid))) stdresid[!is.finite(stdresid)]=0
	# p=1 normal case, p=2 squared std. residuals
	# Q-Statistics on Standardized Residuals
	#H0 : No serial correlation ==> Accept H0 when prob. is High [Q < Chisq(lag)]
	box10 = Weighted.Box.test(stdresid, lag = 1, type = "Ljung-Box", fitdf = 0, if(p==2) sqrd.res = TRUE else sqrd.res = FALSE)
	box15 = Weighted.Box.test(stdresid, lag = max(2, 2*df+df-1), type = "Ljung-Box", fitdf = df, if(p==2) sqrd.res = TRUE else sqrd.res = FALSE)
	box20 = Weighted.Box.test(stdresid, lag = max(5, 4*df+df-1), type = "Ljung-Box", fitdf = df, if(p==2) sqrd.res = TRUE else sqrd.res = FALSE)
	LBSR<-matrix(NA,ncol=2,nrow=3)
	LBSR[1:3,1] = c(box10$statistic[[1]],box15$statistic[[1]],box20$statistic[[1]])
	LBSR[1:3,2] = c(box10$p.value[[1]],box15$p.value[[1]],box20$p.value[[1]])
	rownames(LBSR) = c(paste("Lag[1]",sep=""), paste("Lag[2*(p+q)+(p+q)-1][",max(2, 2*df+df-1),"]",sep=""), paste("Lag[4*(p+q)+(p+q)-1][",max(5, 4*df+df-1),"]",sep=""))
	colnames(LBSR) = c("statistic","p-value")
	return(LBSR)
}


.archlmtest = function (x, lags, demean = FALSE)
{
	if(any(!is.finite(x))) x[!is.finite(x)] = 0
	x = as.vector(x)
	if(demean) x = scale(x, center = TRUE, scale = FALSE)
	lags = lags + 1
	mat = embed(x^2, lags)
	arch.lm = summary(lm(mat[, 1] ~ mat[, -1]))
	STATISTIC = arch.lm$r.squared * length(resid(arch.lm))
	names(STATISTIC) = "Chi-squared"
	PARAMETER = lags - 1
	names(PARAMETER) = "df"
	PVAL = 1 - pchisq(STATISTIC, df = PARAMETER)
	METHOD = "ARCH LM-test"
	result = list(statistic = STATISTIC, parameter = PARAMETER,
			p.value = PVAL, method = METHOD)
	class(result) = "htest"
	return(result)
}


.weightedarchlmtest = function (x, sigma, lags, fitdf = 2, demean = FALSE)
{
	if(any(!is.finite(x))) x[!is.finite(x)] = 0
	x = as.vector(x)
	if(demean) x = scale(x, center = TRUE, scale = FALSE)
	result = Weighted.LM.test(x, sigma^2, lag = lags, type = c("correlation", "partial")[1], fitdf = fitdf, weighted=TRUE) 
	return(result)
}


.nyblomTest = function(object)
{
	#pnames = rownames(object@fit$ipars[object@fit$ipars[,"Estimate"]==1,])
	grad = object@fit$scores
	pnames = colnames(grad)
	if(is.null(pnames)) pnames = rownames(object@fit$ipars[object@fit$ipars[,"Estimate"]==1,,drop=FALSE])
	# special case when fixed parameters exist but the fixed.se was used in fit.control options
	#if( length(pnames) != dim(grad)[2] && dim(grad)[2] == length(object@fit$ipars[object@fit$ipars[,"Include"]==1,]) ){
	#	pnames = rownames(object@fit$ipars[object@fit$ipars[,"Include"]==1,])
	#}
	if(is(object, "ARFIMAfit")) res = object@fit$residuals/object@model$pars["sigma", 1] else res = object@fit$residuals/object@fit$sigma
	res[!is.finite(res) | is.nan(res) | is.na(res)] = 0
	nn = length(res)
	hes = t(grad)%*%(grad)
	shes = try(solve(hes), silent = TRUE)
	if(inherits(shes,"try-error")){
		shes = try( solve(object@fit$hessian) )
		if(inherits(shes,"try-error")){
			IndividualStat = matrix(rep(NA, length(pnames)), ncol=1)
			IndividualCritical = .nyblomCritical(1)
			JointCritical = .nyblomCritical(length(pnames))
			rownames(IndividualStat) = pnames
			JointStat = NA
		} else{
			zx = matrix(cbind(res, res^2, res^3, res^4), ncol = 4)
			zs = apply(zx, 2, FUN = function(x) scale(x, center = TRUE, scale = FALSE))
			x = as.matrix(apply(grad, 2, FUN = function(x) cumsum(x)))
			xx = t(x)%*%x
			nyblomj = sum(diag(xx%*%shes))/nn
			nyblomt = diag(xx)/(diag(hes)*nn)
			IndividualStat = matrix(nyblomt, ncol=1)
			IndividualCritical = .nyblomCritical(1)
			JointCritical = .nyblomCritical(length(pnames))
			rownames(IndividualStat) = pnames
			JointStat = nyblomj
		}
	} else{
		zx = matrix(cbind(res, res^2, res^3, res^4), ncol = 4)
		zs = apply(zx, 2, FUN = function(x) scale(x, center = TRUE, scale = FALSE))
		x = as.matrix(apply(grad, 2, FUN = function(x) cumsum(x)))
		xx = t(x)%*%x
		nyblomj = sum(diag(xx%*%shes))/nn
		nyblomt = diag(xx)/(diag(hes)*nn)
		IndividualStat = matrix(nyblomt, ncol=1)
		IndividualCritical = .nyblomCritical(1)
		JointCritical = .nyblomCritical(length(pnames))
		rownames(IndividualStat) = pnames
		JointStat = nyblomj
	}
	return(list(IndividualStat=IndividualStat,JointStat=JointStat,
					IndividualCritical=IndividualCritical,JointCritical=JointCritical))
}
.nyblomCritical = function(n){
	# Test for Constancy of Parameters: Hansen Tests (1990 mimeo paper)
	# Null Hypothesis: constant parameters
	cval = matrix(c(0.353, 0.470, 0.748,
					0.610, 0.749, 1.07,
					0.846, 1.01, 1.35,
					1.07,  1.24, 1.60,
					1.28,  1.47, 1.88,
					1.49,  1.68, 2.12,
					1.69,  1.90, 2.35,
					1.89,  2.11, 2.59,
					2.10,  2.32, 2.82,
					2.29,  2.54, 3.05,
					2.49,  2.75, 3.27,
					2.69,  2.96, 3.51,
					2.89,  3.15, 3.69,
					3.08,  3.34, 3.90,
					3.26,  3.54, 4.07,
					3.46,  3.75, 4.30,
					3.64,  3.95, 4.51,
					3.83,  4.14, 4.73,
					4.03,  4.33, 4.92,
					4.22,  4.52, 5.13), byrow=T, ncol=3)
	colnames(cval)=c("10%","5%","1%")
	if(n<=20){
		ans = cval[n,]
		names(ans)=c("10%","5%","1%")
	} else {
		ans="too many parameters"
	}
	return(ans)
}

.signbiasTest = function(object)
{
	if(is(object, "uGARCHfilter")) z = object@filter$z else z = z = object@fit$z
	res = as.numeric(residuals(object))
	z2 = z^2
	n = length(z)
	zminus = as.integer(res<0)
	zplus = 1-zminus
	czminus = zminus*res
	czplus = zplus*res
	cz = cbind(rep(1, n), zminus, czminus, czplus)
	cz = cz[1:(n-1),]
	z2 = matrix(z2[2:n],ncol = 1)
	cm = data.frame(y = z2, const = cz[,1], zminus = cz[,2], czminus = cz[,3], czplus = cz[,4])
	fitA = lm(y~const+zminus+czminus+czplus-1, data = cm)
	resA = residuals(fitA)
	sgmatrix = matrix(ncol = 3, nrow = 4)
	rownames(sgmatrix) = c("Sign Bias","Negative Sign Bias","Positive Sign Bias","Joint Effect")
	colnames(sgmatrix) = c("t-value","prob","sig")
	sgmatrix[1:3,1:2] = abs(summary(fitA)$coef[2:4,3:4])
	jeffect = .linear.hypothesis(fitA, c("zminus = 0", "czminus = 0","czplus  =  0"), test  =  "Chisq")
	sgmatrix[4,1] = jeffect[5]$Chisq[2]
	sgmatrix[4,2] = jeffect[6][2,]
	sgmatrix = as.data.frame(sgmatrix)
	sgmatrix[,3] = .stars(sgmatrix[,2])
	return(sgmatrix)
}

.gofTest = function(object, groups)
{
	modelinc = object@model$modelinc
	idx = object@model$pidx
	if(is(object, "uGARCHfilter")){
		ipars = object@filter$ipars
		z = object@filter$z
	} else{
		ipars = object@fit$ipars
		z = object@fit$z
	}
	# must remove fixed parameters
	dist = object@model$modeldesc$distribution
	cdfv = pdist(dist, q = sort(z), lambda = ipars[idx["ghlambda",1],1], skew = ipars[idx["skew",1],1], 
			shape = ipars[idx["shape",1],1])
	j = length(groups)
	gofmat = matrix(NA, ncol = 3, nrow = j)
	gofmat[,1] = groups
	for(i in 1:j){
		sq = seq(1/groups[i], 1, by = 1/groups[i])
		ni = tabulate(findInterval(cdfv, c(0, sq), rightmost.closed = TRUE, 
						all.inside = FALSE)) 
		ExpValue = length(cdfv)/groups[i]
		gofmat[i,2] = sum( ( ( ni - ExpValue ) ^2 )/ExpValue )
		gofmat[i,3] = pchisq(q=gofmat[i,2], df = groups[i]-1, lower.tail = FALSE)
	}
	colnames(gofmat) = c("group", "statistic", "p-value(g-1)")
	rownames(gofmat) = 1:j
	return(gofmat)
}

# returns various loss functions (vectors)
lossfn.var = function(realized, forecast)
{
	x = cbind(realized, forecast)
	#exc = which(is.na(x))
	#if(length(exc)>0) x = x[-exc,]
	mse  = (x[,1] - x[,2])^2
	# qlike: quasi likelihood (implied by gaussian likelihood)
	qlike = log(x[,2]) + x[,1]/x[,2]
	# r2log: penalizes forecasts asymmetrically in low and high volatility periods
	r2log = ( log(x[,1]/x[,2]) )^2
	mad = abs(x[,1] - x[,2])
	#hmse = (x[,1]*(1/x[,2]))^2
	lossdf = data.frame(cbind(mse, mad, qlike, r2log))
	return(lossdf)
}

lossfn.ret = function(realized, forecast)
{
	x = cbind(realized, forecast)
	#exc = which(is.na(x), arr.ind = TRUE)
	#if(length(exc)>0) x = x[-exc[,1],]
	mse  = (x[,1] - x[,2])^2
	mae = abs(x[,1] - x[,2])
	#hmse = (x[,1]*(1/x[,2]))^2
	# we count zero as positive for investment purposes
	dac = as.integer(.signpluszero(x[,1])==.signpluszero(x[,2]))
	lossdf = data.frame(cbind(mse, mae, dac))
	return(lossdf)
}


.meanlossfn = function(lossdf)
{
	n = dim(lossdf)[1]
	nm = colnames(lossdf)
	ans = apply(lossdf, 2, FUN = function(x) 1/(n-1) * sum(x[!is.na(x)]))
	names(ans) = nm
	return(ans)
}

.medianlossfn = function(lossdf)
{
	n = dim(lossdf)[1]
	xn = dim(lossdf)[2]
	nm = colnames(lossdf)
	ans = apply(lossdf[,-xn], 2, FUN = function(x) median(x, na.rm = TRUE))
	# obviously we cannot have the median of a 0/1 series.
	ans = c(ans, .meanlossfn(lossdf[,xn]))
	names(ans) = nm
	return(ans)
}

# Diebold Mariano Test (JBES 1995)
#.DM.test = function(x, lagq)
#{
	# x = vector of differences btw the loss function of the benchmark model and the loss fn of the competing model
	# lagq = lag considered to calculate the newey-west var-cov matrix
	#n.ahead = forecast horizon
	# Returns:
	# an anova class object from the linear.hypothesis test of library car
	# x = rnorm(100, 0, 0.01) - rnorm(100,0.0001,0.01)
#	exc = which(is.na(x))
#	if(length(exc)>0) x = x[-exc]
#	dm = lm(x~1)
#	dm.test = .linear.hypothesis(dm,  c("(Intercept) = 0"), test = "Chisq", vcov  =  NeweyWest(dm, lag = lagq, 
#					order.by = NULL, prewhite = TRUE, adjust = TRUE))
#	return(dm.test)
#}

# From his 2001 paper. Also added the Jarque Bera Test as reccomended by Dowd
# since the test does not really account for residuals being from the Normal distribution.

BerkowitzTest = function(data, lags = 1, significance = 0.05, tail.test = FALSE, alpha = 0.05)
{
	if(tail.test){
		ans = .BerkowitztLRtail(data, alpha = alpha, significance = significance)
		ans$rho =  NA
		ans$JB = NA
		ans$JBp = NA
		return(ans)
	} else{
		if( lags < 1 ) stop("\nlags must be 1 or greater!") else lags = as.integer(lags)
		x = as.numeric(data) - mean(data)
		n = length(data)
		xlag = NULL
		for(i in 1:lags) xlag = cbind(xlag, .lagx(x, n.lag = i, pad = 0))
		ans = lm(x~xlag-1)
		uLL = sum(dnorm(residuals(ans)[-c(1:lags)], sd = summary(ans)$sigma, log = TRUE))
		rLL = sum(dnorm(data[-c(1:lags)], log = TRUE))
		LR = 2*(uLL - rLL)
		chid = 1-pchisq(LR, 2 + lags)
		if(chid < significance) res = paste("reject NULL") else res = paste("fail to reject NULL")
		H0 = paste("Normal(0,1) with no autocorrelation")
		m1 = sum(x)/n
		xm = (x - m1)
		m2 = sum(xm^2)/n
		m3 = sum(xm^3)/n
		m4 = sum(xm^4)/n
		k1 = (m3/m2^(3/2))^2
		k2 = (m4/m2^2)
		JB = n * k1/6 + n * (k2 - 3)^2/24
		JBp = 1 - pchisq(JB, df = 2)
		return(list(uLL = uLL, rLL = rLL, LR = LR, LRp = chid, H0 = H0, Decision = res, 
						mu = mean(data), sigma = summary(ans)$sigma, rho =  coef(ans)[1:(lags)],
				JB = JB, JBp = JBp))
	}
}

.BerkowitztLRtail = function(data, alpha = 0.05, significance = 0.05){
	
	.lrh = function(pars, x){
		p1 = x[which(x<qnorm(alpha))]
		p2 = x[which(x>=qnorm(alpha))]*0 + qnorm(alpha)
		-( sum(log(dnorm((p1-pars[1])/pars[2])/pars[2])) + sum(log(1-pnorm((p2 - pars[1])/pars[2]) )) )
	}
	tst = solnp(c(0.0, 0.95), fun = .lrh, LB = c(-10, 0.01), UB = c(10, 3), x = data, control=list(trace = FALSE))
	uLL = -tail(tst$values, 1)
	rLL = -.lrh(c(0, 1), data)
	LR = 2 * (uLL - rLL)
	chid = 1 - pchisq(LR, 2)
	if (chid < significance) 
		res = paste("reject NULL")
	else res = paste("fail to reject NULL")
	H0 = paste("Normal(0,1)")
	return(list(uLL = uLL, rLL = rLL, LR = LR, LRp = chid, H0 = H0, 
					Decision = res, mu = tst$par[1], sigma = tst$par[2]))
}

# Tests of Directional Accuracy
DACTest = function(forecast, actual, test = c("PT", "AG"), conf.level = 0.95)
{
  n = length(actual)
  if( length(forecast) != n ) stop("Length of forecast and actual must be the same")
  if( test == "PT"){
    x_t = z_t = y_t = rep(0, n)
    x_t[which(actual>0)] = 1
    y_t[which(forecast>0)] = 1
    p_y = mean(y_t)
    p_x = mean(x_t)
    z_t[which( (forecast*actual)>0 )] = 1
    p_hat = mean(z_t)
    p_star = p_y*p_x + (1 - p_y)*(1-p_x)
    p_hat_var = (p_star*(1-p_star))/n
    p_star_var = ((2*p_y-1)^2*(p_x*(1-p_x)))/n + ((2*p_x-1)^2*(p_y*(1-p_y)))/n + (4*p_x*p_y*(1-p_x)*(1-p_y))/n^2
    s_n = (p_hat - p_star)/sqrt(p_hat_var - p_star_var)
    ans = list(Test = "Pesaran and Timmermann", Stat = s_n, p.value = 1-pnorm(s_n), H0 = "Independently Distributed",
    Decision = if( s_n >  qnorm(conf.level) ) "Reject  H0" else "Fail to Reject H0",
    DirAcc = p_hat)
  } else{
    r_t=sign(forecast)*actual
    A_t = mean(r_t)
    B_t = mean(sign(forecast))*mean(actual)
    p_y = 0.5 * (1 + mean(sign(forecast)))
    
    V_EP = (4/(n^2))*p_y*(1-p_y)*sum((actual-mean(actual))^2)
    EP = (A_t-B_t)/sqrt(V_EP)
    ans = list(Test = "Anatolyev and Gerko", Stat = EP, p.value = 1-pnorm(EP), H0 = "No Predictability",
    Decision = if( EP >  qnorm(conf.level) ) "Reject  H0" else "Fail to Reject H0",
    DirAcc = sum(r_t>0)/n)
  }
  return( ans )
}

VaRTest = function(alpha = 0.05, actual, VaR, conf.level = 0.95){
	N = length(actual)
	VaRn = floor(N * alpha)
	if(N != length(VaR)) stop("\nlength of realized not equal to length of VaR!")
	tmp = LR.cc.test(p = alpha, actual = actual, VaR = VaR, conf.level = conf.level)
	ans = list()
	ans$expected.exceed = floor(alpha*tmp$TN)
	ans$actual.exceed = tmp$N
	ans$uc.H0 = "Correct Exceedances"
	ans$uc.LRstat = tmp$stat.uc
	ans$uc.critical = tmp$crit.val.uc
	ans$uc.LRp = tmp$p.value.uc
	ans$uc.Decision = ifelse(ans$uc.LRp<(1-conf.level), "Reject H0", "Fail to Reject H0")
	
	ans$cc.H0 = "Correct Exceedances & Independent"
	ans$cc.LRstat = tmp$stat.cc
	ans$cc.critical = tmp$crit.val.cc
	ans$cc.LRp = tmp$p.value.cc
	ans$cc.Decision = ifelse(ans$cc.LRp<(1-conf.level), "Reject H0", "Fail to Reject H0")
	return(ans)
}

ESTest = function(alpha = 0.05, actual, ES, VaR, conf.level = 0.95, boot = FALSE, n.boot = 1000){
	N = length(actual)
	if(N != length(VaR)) stop("\nlength of realized not equal to length of VaR!")
	idx = which(as.numeric(VaR)>actual)
	z = (ES[idx] - actual[idx])
	.fn = function(x){
		n = length(x)
		1-pnorm(mean(x)/((sqrt(sd(x)^2*((n-1)/n)))/sqrt(n-1)))
	}
	if(boot){
		rbt = matrix(sample(z, size = length(z) * n.boot, replace = TRUE), nrow = n.boot)
		pv = mean( apply(rbt, 1, FUN = function(x) .fn(x)) )
		npv = .fn(z)
	} else{
		npv = .fn(z)
		pv = NA
	}
	ans = list()
	ans$expected.exceed = floor(alpha*N)
	ans$actual.exceed = length(idx)
	# conditional expected shortfall is systematically underestimated
	ans$H1 = "Mean of Excess Violations of VaR is greater than zero"
	ans$boot.p.value = pv
	ans$p.value = npv
	ans$Decision = ifelse(npv<(1-conf.level),"Reject H0", "Fail to Reject H0")
	return(ans)
}


VaRDurTest = function(alpha, actual, VaR, conf.level = 0.95){
	VaR.ind = ifelse(actual < VaR, 1, 0)
	N = sum(VaR.ind)
	TN = length(VaR.ind)
	D = diff(which(VaR.ind==1))
	C = rep(0, length(D))
	# left-censored
	if(VaR.ind[1]==0){
		C = c(1, C)
		# the number of days until we get the first hit
		D = c(which(VaR.ind==1)[1], D)
	}
	# right-censored
	if(VaR.ind[TN]==0){
		C = c(C, 1)
		# the number of days after the last one in the hit sequence
		D = c(D, TN - tail(which(VaR.ind==1), 1))
	}
	N = length(D)
	sol = try(optim(par = 2, fn = .likDurationW, gr = NULL, D = D, C = C, N = N, 
			method = "L-BFGS-B", lower = 0.001, upper = 10, control = list(trace=0)), silent = TRUE)
	b = sol$par
	uLL = -sol$value
	rLL = -.likDurationW(1, D, C, N)
	LR = 2*(uLL - rLL)
	LRp = 1 - pchisq(LR, 1)
	H0 = "Duration Between Exceedances have no memory (Weibull b=1 = Exponential)"
	#i.e. whether we fail to reject the alternative in the LR test that b=1 (hence correct model)
	Decision = ifelse(LRp<(1-conf.level),"Reject H0", "Fail to Reject H0")
	return(list(b = b, uLL = uLL, rLL = rLL, LRp = LRp, H0 = H0, Decision = Decision))
}

.likDurationW = function(pars, D, C, N){
	b = pars[1]
	a = ( (N - C[1] - C[N])/(sum(D^b)) )^(1/b)
	lik = C[1]*log(.pweibull(D[1],a,b,survival=TRUE)) + (1-C[1])*.dweibull(D[1], a, b, log = TRUE) + 
					sum(.dweibull(D[2:(N-1)], a, b, log = TRUE) ) + C[N]*log(.pweibull(D[N],a,b,survival=TRUE) )  + 
									(1 - C[N]) *.dweibull(D[N], a, b, log = TRUE)
	if(!is.finite(lik) || is.nan(lik)) lik = 1e10 else lik = -lik
	return(lik)
}
	
# When b=1 we get the exponential
.dweibull = function(D, a, b, log = FALSE){
	# density of Weibull
	pdf = b * log(a) + log(b) + (b - 1) * log(D) - (a * D)^b
	if(!log) pdf = exp(pdf)
	return(pdf)
}
.pweibull = function(D, a, b, survival = FALSE){
	# distribution of Weibull
	cdf = 1 - exp(-(a*D)^b)
	if(survival) cdf = 1 - cdf
	return(cdf)
}
.hweibull = function(D, a, b){
	# hazard of Weibull
	h = (a^b)*b*(D^(b-1))
	return(h)
}

#####################################################################################
# CHANGELOG (15-08-2012): Changed kurt argument in GMMTest to default of 3 (Normal) 
# ...was badly set at 0. 

GMMTest = function(z, lags = 1, skew=0, kurt=3, conf.level = 0.95){
	if(length(skew)>1) sk = skew[-c(1:lags)] else sk = skew
	if(length(kurt)>1) ku = kurt[-c(1:lags)] else ku = kurt
	z = matrix(z, ncol = 1)
	N = dim(z)[1] - lags
	zlag = z[-c(1:lags), , drop = FALSE]
	orthmat = matrix(NA, ncol = 8, nrow = 3)
	colnames(orthmat) = c("E[z]", "E[z^2]-1", "E[z^3]", "E[z^4]-3", "Q2", "Q3", "Q4","J")
	rownames(orthmat) = c("mean", "var", "t.value")
	f1 = zlag[,1]
	orthmat[1:3,1] = c(mean(f1), mean(f1^2)/N, mean(f1)/sqrt(mean(f1^2)/N))
	f2 = (zlag[,1]^2)-1
	orthmat[1:3,2] = c(mean(f2), mean(f2^2)/N, mean(f2)/sqrt(mean(f2^2)/N))
	f3 = zlag[,1]^3-sk
	orthmat[1:3,3] = c(mean(f3), mean(f3^2)/N, mean(f3)/sqrt(mean(f3^2)/N))
	f4 = zlag[,1]^4-ku
	orthmat[1:3,4] = c(mean(f4), mean(f4^2)/N, mean(f4)/sqrt(mean(f4^2)/N))
	M = rbind(t(f1), t(f2), t(f3), t(f4))
	
	tmp1 = .waldcomomtest(z[,1]^2-1, lags, N)
	orthmat[3, 5] = tmp1$tval[lags+1]
	
	tmp2 = .waldcomomtest(z[,1]^3-skew, lags, N)
	orthmat[3, 6] = tmp2$tval[lags+1]
	
	tmp3 = .waldcomomtest(z[,1]^4-kurt, lags, N)
	orthmat[3, 7] = tmp3$tval[lags+1]
	
	M = rbind(M, tmp1$h, tmp2$h, tmp3$h)
	g = c(as.numeric(orthmat[1,1:4]), tmp1$g, tmp2$g, tmp3$g)
	
	# all moments
	S = (M %*% t(M))/N
	jtval = N*t(g)%*%solve(S)%*%g
	orthmat[3,8] = jtval
	p = rep(conf.level, 2)
	df = c(lags, lags, lags, (4+3*lags))	
	critical.values = qchisq(p,df)
	# i.e. if tval>critical value reject the NULL
	Decision = NULL
	Decision[1] = ifelse(orthmat[3, 5]<critical.values[1], "Fail to Reject H0", "Reject H0")
	Decision[2] = ifelse(orthmat[3, 6]<critical.values[2], "Fail to Reject H0", "Reject H0")
	Decision[3] = ifelse(orthmat[3, 7]<critical.values[3], "Fail to Reject H0", "Reject H0")
	Decision[4] = ifelse(orthmat[3, 8]<critical.values[4], "Fail to Reject H0", "Reject H0")
	
	H0 = "[Q-Moment Conditions] Model is Correctly Specified"
	moment.mat = orthmat[1:3,1:4]
	joint.mat = rbind(orthmat[3,5:8], critical.values)
	rownames(joint.mat) = c("t-value", "critical.value")
	return(list(joint.mat = joint.mat, moment.mat = moment.mat, H0 = H0, Decision = Decision) )
}

.waldcomomtest = function(z, lags, N){
	f0 = z %*% matrix(1, nrow = 1, ncol = lags)
	fx = NULL
	for(i in 1:lags) fx = cbind(fx, .lagx(z, n.lag = i, pad = 0))
	fx = fx[-c(1:lags), , drop = FALSE]
	f0 = f0[-c(1:lags), , drop = FALSE]
	fmat = f0 * fx
	g = colMeans(fmat)
	varg = apply(fmat, 2, FUN = function(x) mean(x^2)/N)
	tval = g/varg
	h = t(fmat)
	S = (h%*%t(h))/N
	joint = N*t(g)%*%solve(S)%*%g
	tval = c(tval, joint)
	return(list(tval = tval, h = h, g = g, varg = varg, S = S))
}

# currenly only quartic kernel implemented
# Hong and Li Test
# M(1,2) test for ARCH-in-Mean 
# M(2,1) tests for leverage
HLTest = function(PIT, lags = 4, kernel = "quartic", conf.level = 0.95){
	p = lags
	Vcon = 2*(50/49 - 300/294 + 1950/1960 - 900/1568 + 450/2304)^2
	Qhatvector = matrix(0, p,1)
	res = rep(0, 7)
	Acon2 = integrate(funb, 0, 1)$value
	T = length(PIT)
	hpit = sd(PIT) * T^(-1/6)
	Acon11 = (1/hpit-2)*(5/7)
	Acon_1 = (Acon11 + 2*Acon2)^2 - 1
	for(i in 1:p){
		Mhat1 = gauss_legendre2D(f = ghat, 0,1,0,1, pit = PIT, i = i, T = T, hpit = hpit) 			
		Qhatvector[i,1] = ( (T-i)*hpit*Mhat1 - hpit*Acon_1 )/sqrt(Vcon)
	}
	res[7] = sum(Qhatvector[,1])/sqrt(p)
	
	res[1] = Momstat(1,1, p, PIT, T)
	res[2] = Momstat(2,2, p, PIT, T)
	res[3] = Momstat(3,3, p, PIT, T)
	res[4] = Momstat(4,4, p, PIT, T)
	res[5] = Momstat(1,2, p, PIT, T)
	res[6] = Momstat(2,1, p, PIT, T)
	
	names(res) = c("M(1,1)", "M(2,2)", "M(3,3)", "M(4,4)", "M(1,2)", "M(2,1)", "W")
	Decision = NULL
	RejectH0 = as.logical( res>rep( qnorm(conf.level), 7 ))
	Decision[1] = ifelse(RejectH0[1], "Reject H0", "Fail to Reject H0")
	Decision[2] = ifelse(RejectH0[2], "Reject H0", "Fail to Reject H0")
	Decision[3] = ifelse(RejectH0[3], "Reject H0", "Fail to Reject H0")
	Decision[4] = ifelse(RejectH0[4], "Reject H0", "Fail to Reject H0")
	Decision[5] = ifelse(RejectH0[5], "Reject H0", "Fail to Reject H0")
	Decision[6] = ifelse(RejectH0[6], "Reject H0", "Fail to Reject H0")
	Decision[7] = ifelse(RejectH0[7], "Reject H0", "Fail to Reject H0")
	names(Decision) = c("M(1,1)", "M(2,2)", "M(3,3)", "M(4,4)", "M(1,2)", "M(2,1)", "W")
	H0 = NULL
	H0[1] = "M(1,1) Correctly Specified"
	H0[2] = "M(2,2) Correctly Specified"
	H0[3] = "M(3,3) Correctly Specified"
	H0[4] = "M(4,4) Correctly Specified"
	H0[5] = "M(1,2) Correctly Specified"
	H0[6] = "M(2,1) Correctly Specified"
	H0[7] = "Model Correctly Specified"
	return(list(statistic = res, H0 = H0, Decision = Decision))
}

funb = function(b){
	tmp1 = (8/15 + b - (2/3)*(b^3)+(1/5)*(b^5))^(-2)
	tmp2 = b*( (1-b^2)^4 )+128/315+(8/3)*(b^3)-(24/5)*(b^5)+(24/7)*(b^7)-(8/9)*(b^9)
	ans = tmp1 * tmp2
	return(ans)
}

ghat = function(x, pit, i, T, hpit){
	z1 = x[1]
	z2 = x[2]
	n = length(pit)
	b1 = fnbound( z1, pit[-c(1:i)], hpit)
	b2 = fnbound( z2, pit[-c((n-i+1):n)], hpit)
	g = sum(b1* b2)/(T-i)
	return( (g - 1)^2 )
}

fnbound = function(x, y, hpit)
{
	K1 = quartic( (x - y)/hpit )/hpit
	if( (x >= 0)*(x < hpit) ){
		K2 = integrate(quartic, (-x/hpit), 1)$value
		return(K1/K2)
	} else if ((x >= hpit)*(x <= (1-hpit))){
		return( K1 )
	} else if( (x>(1-hpit))*(x <= 1) ){
		K2 = integrate(quartic, -1, (1-x)/hpit )$value
		return(K1/K2)
	}
}

quartic = function(z){
	kernell = rep(0, length(z))
	nonzeros = which(z<=1 & z>=-1)
	v = z[nonzeros]
	v = 1-v^2
	kernell[nonzeros] = (15/16)*(v^2)
	return(kernell)
}


Momstat = function(m, ll, p, pit, T){
	part1 = 0
	part2 = 0
	part3 = 0
	for(j in 1:p){
		w = 1 - j/p
		part1 = part1 + (w^2)*(T-j)*((Rcc1(j, m, ll, pit, T)/Rcc1(0, m, ll, pit, T))^2)
		part2 = part2 + w^2
		part3 = part3 + w^4
	}
	return( (part1-part2)/part3 )
}

Rcc1 = function(j, m, ll, pit, T){
	R1 = sum( (pit[(j+1):T]^m) * (pit[1:(T-j)]^ll) ) / T
	R2 = sum( (pit[(j+1):T]^m) ) / T
	R3 = sum( (pit[1:(T-j)]^ll) ) / T
	return( R1 - R2*R3 )
}


##############################################################################
# 2D Integration
# Code from forums:
# http://tolstoy.newcastle.edu.au/R/e5/help/08/09/2919.html
# posted by Earl F. Glynn

gauss_legendre2D_helper <- function(f, x, a2,b2, nodes, weights, ...) {
	C <- (b2 - a2) / 2
	D <- (b2 + a2) / 2
	sum <- 0.0
	for (i in 1:length(nodes))
	{
		y <- nodes[i]*C + D
		sum <- sum + weights[i] * f(c(x,y), ...)
	}
	return(C * sum)
}

gauss_legendre2D <- function(f, a1,b1, a2,b2, ...) {
	#library(statmod)
	#N <- 12
	#GL <- gauss.quad(N)
	nodes = c(-0.9815606, -0.9041173, -0.7699027, -0.5873180, -0.3678315, 
			-0.1252334, 0.1252334, 0.3678315, 0.5873180, 0.7699027, 0.9041173, 
			0.9815606)
	weights = c(0.04717534, 0.10693933, 0.16007833, 0.20316743, 
			0.23349254, 0.24914705, 0.24914705, 0.23349254, 0.20316743, 
			0.16007833, 0.10693933, 0.04717534)	
	C <- (b1 - a1) / 2
	D <- (b1 + a1) / 2
	sum <- 0.0
	for (i in 1:length(nodes))
	{
		x <- nodes[i]*C + D
		sum <- sum + weights[i] * gauss_legendre2D_helper(f, x, a2, b2, nodes, weights, ...)
	}
	
	return(C * sum)
}
#####################################################################################



.meanlossfn = function(lossdf)
{
	n = dim(lossdf)[1]
	nm = colnames(lossdf)
	ans = apply(lossdf, 2, FUN = function(x) 1/(n-1) * sum(x[!is.na(x)]))
	names(ans) = nm
	return(ans)
}

.medianlossfn = function(lossdf)
{
	n = dim(lossdf)[1]
	xn = dim(lossdf)[2]
	nm = colnames(lossdf)
	ans = apply(lossdf[,-xn], 2, FUN = function(x) median(x, na.rm = TRUE))
	# obviously we cannot have the median of a 0/1 series.
	ans = c(ans, .meanlossfn(lossdf[,xn]))
	names(ans) = nm
	return(ans)
}



.signpluszero = function(x)
{
	z = as.numeric(x)
	ans = rep(0, length(z))
	ans[z>=0] = 1
	ans
}


.VaRplot = function(varname , p, actual, dates, VaR)
{
	
	y.actual = actual
	y.VaR = VaR
	x.dates = dates
	zd = .makedate(as.character(dates))
	if(zd$status == 0) dates = 1:length(dates) else dates = zd$dates
	title = paste("Daily Returns and Value-at-Risk Exceedances\n","(Series: ", varname,", alpha=", p,")",sep="")
	if(zd$status){
		plot(x = as.Date(dates, format = zd$dformat), y = y.actual, type = "n",
				main = title, ylab = "Daily Log Return", xlab = "time", 
				ylim = c(min(y.actual, y.VaR), max(y.actual, y.VaR)))
	} else{
		plot(as.numeric(dates), y.actual, type = "n",
				main = title, ylab = "Daily Log Return", xlab = "time", 
				ylim = c(min(y.actual, y.VaR), max(y.actual, y.VaR)))
	}
	abline(h = 0, col = 2, lty = 2)
	sel  =  which(y.actual>0)
	points(dates[sel], y.actual[sel], pch = 18, col = "green")
	sel  =  which(y.actual<0)
	points(dates[sel], y.actual[sel], pch = 18, col = "orange")
	sel  =  which(y.actual<y.VaR)
	points(dates[sel], y.actual[sel], pch = 18, col = "red")
	if(zd$status){
		lines(x = as.Date(dates, format = "%Y-%m-%d"), y = y.VaR, lwd = 2, col = "black")
	} else{
		lines(as.numeric(dates), y.VaR, lwd = 2, col = "black")
	}
	legend("topleft", max(actual),c("return >= 0","VaR <= return < 0","return < VaR","VaR"),
			col=c("green","orange","red","black"), cex=0.75,
			pch = c(18,18,18,-1), lty=c(0,0,0,1), lwd=c(0,0,0,2), bty = "n")
}



.VaRreport = function(varname, garchmodel, distribution, p, actual, VaR, conf.level = 0.95)
{
	actual<-as.matrix(actual)
	VaR<-as.matrix(VaR)
	result <- LR.cc.test(p=p,actual=actual,VaR=VaR,conf.level=conf.level)
	cat("VaR Backtest Report\n")
	cat("===========================================\n")
	cat(paste("Model:\t\t\t\t","",garchmodel,"-",distribution,"\n",sep=""))
	cat(paste("Backtest Length:\t",result$TN,"\n",sep=""))
	cat(paste("Data:\t\t\t\t","",varname,"\n",sep=""))
	cat("\n==========================================\n")		
	cat(paste("alpha:\t\t\t\t",round(100*p,1),"%\n",sep=""))
	cat(paste("Expected Exceed:\t",round(p*result$TN,1),"\n",sep=""))
	cat(paste("Actual VaR Exceed:\t",result$N,"\n",sep=""))
	cat(paste("Actual %:\t\t\t",round(100*result$N/result$TN,1),"%\n",sep=""))
	cat("\nUnconditional Coverage (Kupiec)")
	cat("\nNull-Hypothesis:\tCorrect Exceedances\n")
	cat(paste("LR.uc Statistic:\t",round(result$stat.uc,3),"\n",sep=""))
	cat(paste("LR.uc Critical:\t\t",round(result$crit.val.uc,3),"\n",sep=""))
	cat(paste("LR.uc p-value:\t\t",round(result$p.value.uc,3),"\n",sep=""))
	cat(paste("Reject Null:\t\t",ifelse(round(result$p.value.uc,3)< (1-conf.level), "YES","NO"),"\n",sep=""))
	cat("\nConditional Coverage (Christoffersen)")
	cat("\nNull-Hypothesis:\tCorrect Exceedances and\n\t\t\t\t\tIndependence of Failures\n")
	cat(paste("LR.cc Statistic:\t",round(result$stat.cc,3),"\n",sep=""))
	cat(paste("LR.cc Critical:\t\t",round(result$crit.val.cc,3),"\n",sep=""))
	cat(paste("LR.cc p-value:\t\t",round(result$p.value.cc,3),"\n",sep=""))
	cat(paste("Reject Null:\t\t",ifelse(result$reject,"YES","NO"),"\n",sep=""))
}


########################################################################
# Available in a number of locations/textbooks.
# This code originally presented in a webcast by Guy Yollin Feb 2006 (from Insightful)
# Functions to perform Hypothesis test
# on VaR models based on # of exceedances
# calc LR.uc statistic
LR.cc.test = function (p, actual, VaR, conf.level = 0.95) 
{
	result = .LR.cc(p = p, actual = actual, VaR = VaR)
	crit.val.uc = qchisq(conf.level, df = 1)
	crit.val.cc = qchisq(conf.level, df = 2)
	p.value.cc = 1 - pchisq(result$stat.cc, df = 2)
	p.value.uc = 1 - pchisq(result$stat.uc, df = 1)
	reject = ifelse(p.value.cc < 1 - conf.level, TRUE, FALSE)
	return(list(stat.cc = result$stat.cc, stat.uc = result$stat.uc, 
					p.value.cc = p.value.cc, p.value.uc = p.value.uc, conf.level = conf.level, 
					reject = reject, N = result$N, TN = result$TN, crit.val.uc = crit.val.uc, 
					crit.val.cc = crit.val.cc))
}

.LR.cc = function (p, actual, VaR) 
{
	VaR.ind = ifelse(actual < VaR, 1, 0)
	N = sum(VaR.ind)
	TN = length(VaR.ind)
	T00 = sum(c(0, ifelse(VaR.ind[2:TN] == 0 & VaR.ind[1:(TN - 1)] == 0, 1, 0)))
	T11 = sum(c(0, ifelse(VaR.ind[2:TN] == 1 & VaR.ind[1:(TN - 1)] == 1, 1, 0)))
	T01 = sum(c(0, ifelse(VaR.ind[2:TN] == 1 & VaR.ind[1:(TN - 1)] == 0, 1, 0)))
	T10 = sum(c(0, ifelse(VaR.ind[2:TN] == 0 & VaR.ind[1:(TN - 1)] == 1, 1, 0)))
	T0 = T00 + T01
	T1 = T10 + T11
	pi0 = T01/T0
	pi1 = T11/T1
	pe = (T01 + T11)/(T0 + T1)
	# stat.ind = -2 * log((1 - pe)^(T00 + T10) * pe^(T01 + T11)) + 2 * log((1 - pi0)^T00 * pi0^T01 * (1 - pi1)^T10 * pi1^T11)
	stat.ind = -2 *( (T00 + T10)*log(1 - pe) + (T01 + T11)*log(pe)) + 2 * (T00*log(1 - pi0)+T01*log(pi0)+T10*log(1 - pi1)+T11*log(pi1))
	stat.uc = .LR.uc(p = p, TN = TN, N = N)
	stat.cc = stat.uc + stat.ind
	return(list(stat.cc = stat.cc, stat.uc = stat.uc, N = N, 
					TN = TN))
}

.LR.uc = function (p, TN, N) 
{
	stat.uc = -2 *( (TN - N)*log(1 - p)+ N*log(p) ) + 2 * ( (TN - N)*log(1 - N/TN)+N*log(N/TN) )
	return(stat.uc)
}

.Log = function(x){
	ans = log(x)
	#if(!is.finite(ans)) ans = sign(ans) * 1e10
	ans
}

###############################################################################
# MCS Test of Hansen, Lunde and Nason
# Partially transalted from Kevin Sheppard's mcs matlab function
mcsTest = function(losses, alpha, nboot = 100, nblock = 1, boot = c("stationary", "block"))
{
	n = NROW(losses)
	m = NCOL(losses)
	bsdat = bootstrap(1:n, nboot, nblock, type = boot[1])$B
	dijbar = matrix(0, m, m)
	for(j in 1:m) dijbar[j, ] = colMeans(losses - losses[,j])
	dijbarstar = array(0, dim = c(m, m, nboot))
	for(i in 1:nboot){
		tmp = colMeans(losses[bsdat[,i],,drop = FALSE])
		for(j in 1:m){
			dijbarstar[j,,i] = tmp - tmp[j]
		}
	}
	tmp = array(NA, dim = c(m, m, nboot))
	for(i in 1:nboot) tmp[,,i] = dijbar
	xtmp = ((dijbarstar - tmp)^2)
	vardijbar = matrix(0, m, m)
	for(i in 1:nboot) vardijbar = vardijbar + xtmp[,,i]
	vardijbar = vardijbar/nboot
	vardijbar = vardijbar + diag(m)
	tmp2 = array(NA, dim = c(m, m, nboot))
	for(i in 1:nboot) tmp2[,,i] = sqrt(vardijbar)
	z0=(dijbarstar-tmp)/tmp2
	zdata0=dijbar/sqrt(vardijbar)
	excludedR = matrix(0, m,1)
	pvalsR = matrix(1,m,1)
	for(i in 1:(m-1)){
		included = setdiff(1:m,excludedR)
		mx = length(included)
		z = z0[included, included,]
		empdistTR = apply(abs(z), 3, function(x) max(x))
		zdata = zdata0[included,included]
		TR = max(.Call("colMaxRcpp", zdata, PACKAGE = "rugarch"))
		pvalsR[i] = mean(empdistTR>TR)
		dibar = colMeans(dijbar[included,included])*(mx/(mx-1))
		dibstar = apply(dijbarstar[included,included,,drop=FALSE], 3, function(x) colMeans(x))*(mx/(mx-1))
		vardi=colMeans((t(dibstar)-matrix(dibar, nrow = nboot, ncol = mx, byrow = TRUE))^2)
		tx=dibar/sqrt(vardi)
		temp = max(tx)
		modeltoremove = which(tx == max(tx))[1]
		excludedR[i]=included[modeltoremove]
	}
	maxpval=pvalsR[1]
	for(i in 2:m){
		if(pvalsR[i]<maxpval){
			pvalsR[i] = maxpval
		} else{
			maxpval = pvalsR[i]
		}
	}
	excludedR[NROW(excludedR)] = setdiff(1:m,excludedR)
	pl = which(pvalsR >= alpha)[1]
	includedR = excludedR[pl:m]
	if(pl==1) excludedR = NULL else excludedR = excludedR[1:(pl-1)]
	excludedSQ = matrix(0,m,1)
	pvalsSQ = matrix(1,m,1)
	for(i in 1:(m-1)){
		included = setdiff(1:m,excludedSQ)
		mx = length(included)
		z = z0[included,included,]
		empdistTSQ = apply(z^2, 3, function(x) sum(x))/2
		zdata = zdata0[included,included]
		TSQ = sum(colSums(zdata^2))/2
		pvalsSQ[i] = mean(empdistTSQ>TSQ)
		dibar = colMeans(dijbar[included,included])*(mx/(mx-1))
		dibstar = apply(dijbarstar[included,included,,drop=FALSE], 3, function(x) colMeans(x))*(mx/(mx-1))
		vardi = colMeans((t(dibstar)-matrix(dibar, nrow = nboot, ncol = mx, byrow = TRUE))^2)
		tx = dibar/sqrt(vardi)
		temp = max(tx)
		modeltoremove = which(tx == max(tx))[1] 
		excludedSQ[i]=included[modeltoremove]
	}
	maxpval = pvalsSQ[1]
	for( i in 2:m){
		if(pvalsSQ[i]<maxpval){
			pvalsSQ[i] = maxpval
		} else{
			maxpval = pvalsSQ[i]
		}
	}
	excludedSQ[NROW(excludedSQ)] = setdiff(1:m,as.numeric(excludedSQ))
	pl=which(pvalsSQ>=alpha)[1]
	includedSQ = excludedSQ[pl:m]
	if(pl==1) excludedSQ = NULL else excludedSQ = excludedSQ[1:(pl-1)]	
	return(list(includedR = includedR, pvalsR = pvalsR, excludedR = excludedR, 
					includedSQ = includedSQ, pvalsSQ = pvalsSQ, excludedSQ = excludedSQ))
}

bootstrap = function(data, nboot = 10, nblock = 5, type = c("stationary", "block"))
{
	type = match.arg(type[1], c("stationary", "block"))
	ans = switch(type,
			stationary = .stationary_bootstrap(data, nboot, nblock),
			block = .block_bootstrap(data, nboot, nblock))
	return(ans)
}

.stationary_bootstrap = function(data, nboot = 100, nblock = 1){
	n = length(data)
	p = 1/nblock
	idx = matrix(0, n, nboot)
	idx[1, ] = ceiling(n * matrix(runif(1*nboot), 1, nboot))
	sel_idx = matrix(runif(n*nboot), n, nboot)<p
	cn = sum(colSums(sel_idx))
	idx[sel_idx] = ceiling(matrix(runif(1*cn), 1, cn)*n)
	for(i in 2:n){
		idx[i,!sel_idx[i,]] = idx[i-1,!sel_idx[i,]]+1
	}
	data = c(data, data)
	B = apply(idx, 2, FUN = function(x) data[x])
	return(list(B = B, idx = idx))
}

.block_bootstrap = function(data, nboot, nblock){
	n = length(data)
	data = c(data, data[1:(nblock-1)])
	s = ceiling(n/nblock)
	S = floor(matrix(runif(s*nboot), s, nboot)*n)+1
	idx = matrix(0, s*nblock, nboot)
	A = repmat( 0:(nblock-1), 1, nboot)
	iter=1
	for(i in seq(1, n, by = nblock)){
		idx[i:(i+nblock-1),] = repmat(S[iter, ,drop = FALSE], nblock, 1) + A
		iter = iter + 1
	}
	idx = idx[1:n, ,drop = FALSE]
	B = apply(idx, 2, FUN = function(x) data[x])
	return(list(B = B, idx = idx))
}

repmat = function(a,n,m){
	kronecker(matrix(1,n,m),a)
}

##########################################################################################
# Direct Import of Weighted Tests of FISHER and GALLAGHER (WeightedPortTest package)
Weighted.Box.test = function (x, lag = 1, type = c("Box-Pierce", "Ljung-Box", "Monti"), 
		fitdf = 0, sqrd.res = FALSE, log.sqrd.res = FALSE, abs.res = FALSE, 
		weighted = TRUE)
{
	if(lag<(2*fitdf+fitdf-1)) stop("\nLag must be equal to a minimum of 2*fitdf+fitdf-1")
	if(NCOL(x) > 1) stop("\nx is not a vector or univariate time series")
	if(lag < 1) stop("\nLag must be positive")
	if(fitdf < 0) stop("\nFitdf cannot be negative")
	if((sqrd.res && log.sqrd.res) || (sqrd.res && abs.res) || (log.sqrd.res && abs.res)) stop("Only one option of: sqrd.res, log.sqrd.res or abs.res can be selected")
	DNAME <- deparse(substitute(x))
	type <- match.arg(type)
	if (abs.res) {
		x <- abs(x)
	}
	if (sqrd.res || log.sqrd.res) {
		x <- x^2
	}
	if (log.sqrd.res) {
		x <- log(x)
	}
	if (weighted) {
		if (type == "Monti") {
			METHOD <- "Weighted Monti test (Gamma Approximation)"
			cor <- acf(x, lag.max = lag, type = "partial", plot = FALSE, 
					na.action = na.pass)
			obs <- cor$acf[1:lag]
		}
		else {
			cor <- acf(x, lag.max = lag, type = "correlation", 
					plot = FALSE, na.action = na.pass)
			obs <- cor$acf[2:(lag + 1)]
		}
		if (type == "Ljung-Box") {
			METHOD <- "Weighted Ljung-Box test (Gamma Approximation)"
		}
		n <- sum(!is.na(x))
		weights <- (lag - 1:lag + 1)/(lag)
		if (type == "Box-Pierce") {
			METHOD <- "Weighted Box-Pierce test (Gamma Approximation)"
			STATISTIC <- n * sum(weights * obs^2)
		}
		else {
			STATISTIC <- n * (n + 2) * sum(weights * (1/seq.int(n - 1, n - lag) * obs^2))
		}
		if (sqrd.res) {
			fitdf <- 0
			names(STATISTIC) <- "Weighted X-squared on Squared Residuals for detecting nonlinear processes"
		}
		else if (log.sqrd.res) {
			fitdf <- 0
			names(STATISTIC) <- "Weighted X-squared on Log-Squared Residuals for detecting nonlinear processes"
		}
		else if (abs.res) {
			fitdf <- 0
			names(STATISTIC) <- "Weighted X-squared on Absolute valued Residuals for detecting nonlinear processes"
		}
		else {
			names(STATISTIC) <- "Weighted X-squared on Residuals for fitted ARMA process"
		}
		shape <- (3/4) * (lag + 1)^2 * lag/(2 * lag^2 + 3 * lag + 1 - 6 * lag * fitdf)
		scale <- (2/3) * (2*lag^2 + 3*lag+1 - 6 * lag * fitdf)/(lag*(lag + 1))		
		PARAMETER <- c(shape, scale)
		names(PARAMETER) <- c("Shape", "Scale")
		PVAL <- 1 - pgamma(STATISTIC, shape = shape, scale = scale)
		names(PVAL) <- "Approximate p-value"
	}
	else {
		if (type == "Monti") {
			METHOD <- "Monti test"
			cor <- acf(x, lag.max = lag, type = "partial", plot = FALSE, 
					na.action = na.pass)
			obs <- cor$acf[1:lag]
		}
		else {
			cor <- acf(x, lag.max = lag, type = "correlation", 
					plot = FALSE, na.action = na.pass)
			obs <- cor$acf[2:(lag + 1)]
		}
		if (type == "Ljung-Box") {
			METHOD <- "Ljung-Box test"
		}
		n <- sum(!is.na(x))
		if (type == "Box-Pierce") {
			METHOD <- "Box-Pierce test"
			STATISTIC <- n * sum(obs^2)
		}
		else {
			STATISTIC <- n * (n + 2) * sum((1/seq.int(n - 1, 
										n - lag) * obs^2))
		}
		if (sqrd.res) {
			fitdf <- 0
			names(STATISTIC) <- "X-squared on Squared Residuals for detecting nonlinear processes"
		}
		else if (log.sqrd.res) {
			fitdf <- 0
			names(STATISTIC) <- "X-squared on Log-Squared Residuals for detecting nonlinear processes"
		}
		else if (abs.res) {
			fitdf <- 0
			names(STATISTIC) <- "X-squared on Absolute valued Residuals for detecting nonlinear processes"
		}
		else {
			names(STATISTIC) <- "X-squared on Residuals for fitted ARMA process"
		}
		mydf <- lag - fitdf
		PARAMETER <- c(mydf)
		names(PARAMETER) <- c("df")
		PVAL <- 1 - pchisq(STATISTIC, df = mydf)
		names(PVAL) <- "p-value"
	}
	structure(list(statistic = STATISTIC, parameter = PARAMETER, 
					p.value = PVAL, method = METHOD, data.name = DNAME), 
			class = "htest")
}

Weighted.LM.test <- function (x, h.t, lag = 1, type = c("correlation", "partial"), fitdf = 1, weighted=TRUE) 
{
	### Error Checking
	###
	if (NCOL(x) > 1) 
		stop("x is not a vector or univariate time series");
	if (fitdf >= lag)
		stop("Lag must exceed fitted degrees of freedom");
	if (fitdf < 1)
		stop("Fitted degrees of freedom must be positive");
	if( !(length(x)==length(h.t)) )
		stop("Length of x and h.t must match");
	
	DNAME <- deparse(substitute(x))
	type <- match.arg(type)
	
	x <- x^2/h.t
	
	if( type == "partial") {
		cor <- acf(x, lag.max = lag, type="partial", plot=FALSE, na.action=na.pass)
		obs <- cor$acf[1:lag];
	}
	else {
		cor <- acf(x, lag.max = lag, type="correlation", plot=FALSE, na.action=na.pass)
		obs <- cor$acf[2:(lag + 1)];
	}
	
	
	if(type == "correlation" && weighted) {
		METHOD <- "Weighted Li-Mak test on autocorrelations (Gamma Approximation)"
	}
	else if(type == "partial" && weighted) {
		METHOD <- "Weighted Li-Mak test on partial autocorrelations (Gamma Approximation)"
	}
	else if(type == "correlation" && !weighted) {
		METHOD <- "Li-Mak test on autocorrelations (Chi-Squared Approximation)"
	}
	else {
		METHOD <- "Li-Mak test on partial autocorrelations (Chi-Squared Approximation)"
	}
	
	n <- sum(!is.na(x))
	if(weighted) {
		weights <- (lag - (fitdf+1):lag + (fitdf+1) )/lag;
		obs <- obs[(fitdf+1):lag];
		STATISTIC <- n * sum(weights*obs^2);
		names(STATISTIC) <- "Weighted X-squared on Squared Residuals for fitted ARCH process";
		shape <- (3/4)*(lag + fitdf + 1)^2*(lag - fitdf)/(2*lag^2 + 3*lag + 2*lag*fitdf + 2*fitdf^2 + 3*fitdf + 1);
		scale <- (2/3)*(2*lag^2 + 3*lag + 2*lag*fitdf + 2*fitdf^2 + 3*fitdf + 1)/(lag*(lag + fitdf + 1));
		PARAMETER <- c(shape, scale);
		names(PARAMETER) <- c("Shape", "Scale")
	}
	else {
		weights <- rep(1,(lag-fitdf) );
		obs <- obs[(fitdf+1):lag];
		STATISTIC <- n * sum(weights*obs^2);
		names(STATISTIC) <- "X-squared on Squared Residuals for fitted ARCH process"
		shape <- (lag-fitdf)/2;          # Chi-squared df in Gamma form.
		scale <- 2;
		PARAMETER <- c((lag-fitdf));
		names(PARAMETER) <- c("Degrees of Freedom");
	}
	
	PVAL <- 1 - pgamma(STATISTIC, shape=shape, scale=scale)
	names(PVAL) <- "Approximate p-value"
	
	structure(list(statistic = STATISTIC, parameter = PARAMETER, 
					p.value = PVAL, method = METHOD, data.name = DNAME), 
			class = "htest")
}
# End Import
##########################################################################################
