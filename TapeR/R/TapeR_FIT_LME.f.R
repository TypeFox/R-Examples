TapeR_FIT_LME.f <-
function(Id,x,y,knt_x,ord_x,knt_z,ord_z,IdKOVb = "pdSymm", ...){
#   ************************************************************************************************

		## require(nlme)

	    BS_x  	= BSplines(knots = knt_x, ord = ord_x, der = 0, x = x)   ;# fix(BS.x)
	   	BS_z  	= BSplines(knots = knt_z, ord = ord_z, der = 0, x = x)   ;# fix(BS.x)

		X		= BS_x[1:nrow(BS_x),1:ncol(BS_x)-1,drop=F]     ;	# last column is omitted in order to ensure Ey(x=1)=0
	    Z 		= BS_z[1:nrow(BS_z),1:ncol(BS_z)-1,drop=F]     ;

#   ************************************************************************************************
#							LME Modell E[y|X,Z] = X*b_fix + Z*b_rnd + eps
#   ************************************************************************************************

		if(IdKOVb == "pdDiag"){
			fit.lme = lme(y ~ X-1,	random = list(Id = pdDiag(~Z-1)))    #   VAR[b_rnd] from emp KOV[b_rnd]
		}else{
			fit.lme = lme(y ~ X-1,	random = list(Id = pdSymm(~Z-1)),
									control= lmeControl(msMaxIter=500,msMaxEval=500),
									method = "ML" )
		}
#   ************************************************************************************************

		b_fix 		= as.numeric(fixef(fit.lme))
		KOVb_fix    = vcov(fit.lme)

		if(IdKOVb == "pdDiag"){												   # empirische KOVb_rnd
			VARb_rnd = matrix(getVarCov(fit.lme),ncol=ncol(getVarCov(fit.lme)))
			KORb_rnd = cor(ranef(fit.lme))
			KOVb_rnd = sqrt(VARb_rnd)%*%KORb_rnd%*%sqrt(VARb_rnd)
		}else{
			KOVb_rnd = matrix(getVarCov(fit.lme),ncol=ncol(getVarCov(fit.lme)))
			KORb_rnd = cov2cor(KOVb_rnd)
		}

		theta 		= attr(fit.lme$apVar,"Pars")
		KOV_theta 	= matrix(fit.lme$apVar,ncol=length(theta),byrow=T)

		sig2_eps 	= as.numeric(exp(theta["lSigma"])^2)

		par.lme     = list(	knt_x = knt_x, ord_x = ord_x, knt_z = knt_z, ord_z = ord_z,
							b_fix = b_fix, KOVb_fix = KOVb_fix,
							sig2_eps = sig2_eps, dfRes = anova(fit.lme)$denDF,
							KOVb_rnd = KOVb_rnd, theta = theta, KOV_theta = KOV_theta)

#       --------------------------------------------------------------------------------------------
		return(list(fit.lme = fit.lme,par.lme = par.lme))
#       --------------------------------------------------------------------------------------------
	}
