######################################################################################
# Created 18Jan12: This file contains all the posteriors needed for all the 
#	estimation functions.
######################################################################################

######################################################################################
#
#						a Posterior 
#
######################################################################################
#######################################################################################
# log[Posterior(a | mu.a, H00, k)] = a*[log(mu.a)+log(H00)-log(lambda)] - log(Gamma(a+1)) -
#	log(Gamma(a)) + sum_m sum_p [log(Gamma(2ak^m))-log(Gamma(2gamma.mp*ak^m)) -
#	log(Gamma(2(1-gamma.mp)ak^m))+2gamma.mp*ak^m*log(Rmp)+2(1-gamma.mp)*ak^m*log(1-Rmp)]
#######################################################################################
sample.a = function(aval,lambdaval, mu.aval, H00val, kval, mvec, Mval, Rmpvec, gamma.mpvec){

	a=aval
	mat_a=matrix(rep(a,each=2^Mval-1),nrow=length(a),byrow=TRUE)
	
	a_aux01=matrix(2*gamma.mpvec*kval^mvec,ncol=2^Mval-1,nrow=length(a),byrow=TRUE)
	y1=rowSums(lgamma(mat_a*a_aux01))
	
	a_aux02=matrix(2*(1-gamma.mpvec)*kval^mvec,ncol=2^Mval-1,nrow=length(a),byrow=TRUE)
	y2=rowSums(lgamma(mat_a*a_aux02))
	
	a_aux03=matrix(rep(2*kval^mvec,length(a)),nrow=length(a),byrow=TRUE)
	y3= rowSums(lgamma(mat_a*a_aux03))
	
	
	y4=a*(log(mu.aval)+log(H00val)) - lgamma(a+1) - lgamma(a)-a*log(lambdaval) 
	y5=a*sum(2*kval^mvec*gamma.mpvec*log(Rmpvec)+2*kval^mvec*(1-gamma.mpvec)*log(1-Rmpvec))
	
	aPostval = exp(-y1-y2+y3+y4+y5)

	u = runif(1,0,sum(aPostval))
	cum_aPostval = cumsum(aPostval)
	returnVal = length(which(cum_aPostval < u)) + 1

	return(returnVal)
}

######################################################################################
#
#						lambda Posterior 
#
######################################################################################
# Posterior(lambda | mu.l, H00, a) = exp{-lambda/mu.l}*exp{-H00/lambda}/lambda^a
logLambdaPost = function(lambda, mu.lval, H00val, aval){
	
	-lambda/mu.lval - H00val/lambda - aval*log(lambda)
}

######################################################################################
#
#						Rmp Posteriors 
#
######################################################################################

################### Rmp posterior for proportional covariates ####################
# log[Posterior(Rmp | k, H00, a)] = (2gamma.mp*a*k^m-1+sum(delta*RmpIndic))*log(Rmp) + 
#			(2(1-gamma.mp)*a*k^m-1+sum(delta*1_RmpIndic))*log(1-Rmp)-sum(exp(beta.trt*Ind(trt))*H(Ti))
logRmpPost_PHcovs = function(Rmp, RmpFull, H00val, kval, aval, gamma.mpval, betavec, X.mat,
mval, RmpIndic, one_RmpIndic, deltavals, Mvalue, inBinMat, mat01, formula){
	
	RmpFull[formula] = Rmp
	Hvals = calc_H(mat01, H00val, RmpFull, Mvalue, inBinMat)					
	
	(2*gamma.mpval*aval*kval^mval-1+sum(deltavals*RmpIndic))*log(Rmp) + 
	(2*(1-gamma.mpval)*aval*kval^mval-1+sum(deltavals*one_RmpIndic))*log(1-Rmp) - 
	sum(exp(X.mat%*%betavec)*Hvals)
	
}

# log[Posterior(Rmp | k, H00, a)] = (2gamma.mp*a*k^m-1+sum(delta*RmpIndic))*log(Rmp) + 
#			(2(1-gamma.mp)*a*k^m-1+sum(delta*1_RmpIndic))*log(1-Rmp)-sum(exp(beta.trt*Ind(trt))*H(Ti))
logRmpPost_nonPHBA = function(Rmp, RmpFull, H00val, kval, aval, gamma.mpval, betavec, X.mat,
mval, RmpIndic, one_RmpIndic, deltavals, Mvalue, inBinMat, mat01, formula){
	
	RmpFull[formula] = Rmp
	Hvals = calc_H(mat01, H00val, RmpFull, Mvalue, inBinMat)					
	
	(2*gamma.mpval*aval*kval^mval-1+sum(deltavals*RmpIndic))*log(Rmp) + 
	(2*(1-gamma.mpval)*aval*kval^mval-1+sum(deltavals*one_RmpIndic))*log(1-Rmp) - 
	sum(exp(X.mat%*%betavec)*Hvals)
	
}

######################################################################################
#
#						Beta Posteriors 
#
######################################################################################
################### Beta posterior for proportional case ####################
# log[Posterior(beta | mu.beta, sigma.beta, H)] = beta*sum(delta*trtIndicator)-
#						sum[exp{beta*trtIndicator}*H(Ti)] - (beta-mu.beta)^2/2sigma.beta^2
logbetaPost_PH = function(beta, betaFull, deltavals, Xmatrix, Hvals, whichBeta, mu.beta, sigma.beta){
	
	betaFull[whichBeta] = beta
	
	sum(deltavals*(Xmatrix%*%betaFull)) - sum(exp(Xmatrix%*%betaFull)*Hvals) - 
	(beta-mu.beta)^2/(2*sigma.beta^2)
}
# log[Posterior(beta | mu.beta, sigma.beta, H)] = beta*sum(delta*trtIndicator)-
#						sum[exp{beta*trtIndicator}*H(Ti)] - (beta-mu.beta)^2/2sigma.beta^2
logbetaPost_NPHBA = function(beta, betaFull, deltavals, Xmatrix, Hvals, whichBeta, mu.beta, sigma.beta){
	
	betaFull[whichBeta] = beta
	
	sum(deltavals*(Xmatrix%*%betaFull)) - sum(exp(Xmatrix%*%betaFull)*Hvals) - 
	(beta-mu.beta)^2/(2*sigma.beta^2)
}

######################################################################################
#
#						k Posteriors 
#
######################################################################################
logkPost = function(k, gamma.mpval, a, Rmpvals, Rmp.exp, one_Rmp.ex, mvec, mu.k){
	
	k.pow.mvec = k^mvec
	logpost = -k/mu.k + sum(Rmp.exp*log(Rmpvals) + one_Rmp.ex*log(1-Rmpvals) - 
							lgamma(Rmp.exp*k.pow.mvec)+lgamma(one_Rmp.ex*k.pow.mvec)-
							lgamma(Rmp.exp*k.pow.mvec+one_Rmp.ex*k.pow.mvec))
	return(logpost)
}

######################################################################################
#
#						gamma Posteriors 
#
######################################################################################
loggammaPost = function(gammamp, Rmp, kval, aval, cval, dval, mvecval){
	
	k.pow.mvec.timesa = kval^mvecval*aval
	(2*gammamp*k.pow.mvec.timesa-1)*log(Rmp) + (2*(1-gammamp)*k.pow.mvec.timesa-1)*log(1-Rmp) - 
	(cval-1)*log(gammamp) + (dval-1)*log(1-gammamp) - 
	lgamma(2*gammamp*k.pow.mvec.timesa) - lgamma(2*(1-gammamp)*k.pow.mvec.timesa)
}

