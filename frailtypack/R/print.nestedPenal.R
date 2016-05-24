
"print.nestedPenal" <- function (x, digits = max(options()$digits - 4, 6), ...) 
{
	if (!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
		if (x$type == "counting" & x$AG == FALSE) {
			cat("\n      left truncated structure used")
		}
		if (x$AG == TRUE){
			cat("\n      Calendar timescale")
		}
		cat("\n")
	}
	if (!is.null(x$fail)) {
		cat(" frailtyPenal failed.", x$fail, "\n")
		return()
	}
	savedig <- options(digits = digits)
	on.exit(options(savedig))
	coef <- x$coef
	nvar <- length(x$coef)
	
	if (is.null(coef)){
		x$varH<-matrix(x$varH) 
		x$varHIH<-matrix(x$varHIH)
	}
 #AD:     
	if (x$typeof == 0){
		if (x$n.knots.temp < 4){
			cat("\n")
			cat("  The minimum number of knots is 4","\n")
		}
		if (x$n.knots.temp > 20){
			cat("\n")
			cat("  The maximum number of knots is 20","\n")
		}
		
# 		if ((x$indic.Kappa2 == 0) & (x$nst == 1)){
# 			cat(" Kappa2 is not used  \n")
# 		}
	}else{
		if ((x$typeof == 1) & (x$indic.nb.int == 1)) cat("  The maximum number of time intervals is 20","\n")
	}
#AD
	if (x$istop == 1){
		if (!is.null(coef)){ 
			seH <- sqrt(diag(x$varH))[-c(1,2)]
			seHIH <- sqrt(diag(x$varHIH))[-c(1,2)]
			if (x$typeof == 0){
				tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits - 1))
				if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,x$p.global_chisq)
			}else{
				tmp <- cbind(coef, exp(coef), seH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits - 1))
				if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,x$p.global_chisq)
			}	
			cat("\n")
			cat(" Nested Frailty model parameter estimates using ", "\n")

			if (x$typeof == 0){
				cat(" using a Penalized Likelihood on the hazard functions",  "\n")
			}else{
				cat(" using a Parametrical approach for the hazard functions",  "\n")	
			}
			if (x$n.strat>1) cat(" (Stratification structure used) :",x$n.strat,"strata \n")
			
			if (x$typeof == 0){	
				if(x$global_chisq.test==1){
					dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
					
				}
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
				"SE coef (H)", "SE coef (HIH)", "z", "p"))
			}else{
				if(x$global_chisq.test==1){
					dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
					
				}
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
				"SE coef (H)", "z", "p"))
			}
			cat("\n")
			prmatrix(tmp)
			if(x$global_chisq.test==1){
				cat("\n")
				prmatrix(tmpwald)
			}
			cat("\n")
		} 

		alpha <- x$alpha
		temp <- diag(x$varH)[1]
		seH.alpha <- sqrt(((2 * (alpha^0.5))^2) * temp)
		temp <- diag(x$varHIH)[1]
		seHIH.alpha <- sqrt(((2 * (alpha^0.5))^2) * temp)
		
		if (x$noVar1 == 1){
			cat("\n")
			cat("    Nested Frailty model: No covariates \n")
			cat("    -------------------------- \n")
			cat("\n")
		}
		
		cat("  Frailty parameters: \n")
# 		if (x$typeof == 0){
# 			cat("   alpha  (group effect): ", alpha, " (SE(H):", seH.alpha, ")", " (SE(HIH):", seHIH.alpha, ")", sep="", "\n")
# 		}else{
# 			cat("   alpha  (group effect): ", alpha, " (SE(H):", seH.alpha, ")", sep="", "\n")
# 		}
		cat("   alpha  (group effect): ", alpha, " (SE(H):", seH.alpha, ")", "p =", signif(1 - pnorm(alpha/seH.alpha), digits - 1), "\n")

		eta <- x$eta
		temp <- diag(x$varH)[2]
		seH.eta <- sqrt(((2 * (eta^0.5))^2) * temp)
		temp <- diag(x$varHIH)[2]
		seHIH.eta <- sqrt(((2 * (eta^0.5))^2) * temp)
# 		if (x$typeof == 0){
# 			cat("   eta (subgroup effect): ", eta, " (SE(H):", seH.eta, ")", " (SE(HIH):", seHIH.eta, ")", sep="", "\n")
# 		}else{
# 			cat("   eta (subgroup effect): ", eta, " (SE(H):", seH.eta, ")", sep="", "\n")
# 		}
		cat("   eta (subgroup effect): ", eta, " (SE(H):", seH.eta, ")", "p =", signif(1 - pnorm(eta/seH.eta), digits - 1), "\n")

		cat(" \n")

		if (x$typeof == 0){
			cat(paste("    penalized marginal log-likelihood =", round(x$logLikPenal,2)))
			cat("\n")
			cat("    Convergence criteria: \n")
			cat("    parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
			cat("\n")
			cat("    LCV = the approximate likelihood cross-validation criterion\n")
			cat("          in the semi parametrical case     =",x$LCV,"\n")
		}else{
			cat(paste("    marginal log-likelihood =", round(x$logLik,2)))
			cat("\n")
			cat("    Convergence criteria: \n")
			cat("    parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
			cat("\n")
#			cat("    LCV = the approximate likelihood cross-validation criterion\n")
#			cat("          in the parametrical case     =",x$LCV,"\n")
			cat("    AIC = Aikaike information Criterion     =",x$AIC,"\n")
			cat("\n")
			cat("The expression of the Aikaike Criterion is:","\n")
			cat("        'AIC = (1/n)[np - l(.)]'","\n")
			if (x$typeof == 2){
				cat("\n")
				if (x$n.strat == 1){
					cat("      Scale for the weibull hazard function is :",round(x$scale.weib[1],2),"\n")	
					cat("      Shape for the weibull hazard function is :",round(x$shape.weib[1],2),"\n")
				}else{
					cat("      Scale for the weibull hazard function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),"\n")	
					cat("      Shape for the weibull hazard function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),"\n")
				}
					cat("\n")
					cat("The expression of the Weibull hazard function is:","\n")
					cat("        'lambda(t) = (shape.(t^(shape-1)))/(scale^shape)'","\n")
					cat("The expression of the Weibull survival function is:","\n")
					cat("        'S(t) = exp[- (t/scale)^shape]'")
					cat("\n")
			}	
		}
#AD:
		cat("\n")
		cat("    n=", x$n)
		if (length(x$na.action)){ 
			cat("  (", length(x$na.action), " observation deleted due to missing) \n")
		}else{
			cat("\n")
		}
		cat("    n events=", x$n.events, " n groups=", x$groups)
		cat( "\n")
		cat("    number of iterations: ", x$n.iter,"\n")
		if ((x$typeof == 1) & (x$indic.nb.int == 1)){
			cat("      Exact number of time intervals used: 20","\n")
		 }else{
		 	if (x$typeof == 1) cat("    Exact number of time intervals used: ",x$nbintervR,"\n")
		 }
		if (x$typeof == 0){
			cat("    Exact number of knots used: ", x$n.knots, "\n")
			cat("    Value of the smoothing parameter: ", x$kappa, sep=" ")
			cat(", DoF: ", formatC(-x$DoF, format="f",digits=2))
		}
	}else{
		cat("\n")
		cat(" Nested Frailty model parameter estimates using ", "\n")

		if (x$typeof == 0){
			cat(" using a Penalized Likelihood on the hazard functions",  "\n")
		}else{
			cat(" using a Parametrical approach for the hazard functions",  "\n")	
		}
		if (x$n.strat>1) cat(" (Stratification structure used)", "\n") 	
	
		if (x$noVar1 == 1){
			cat("\n")
			cat("    Nested Gamma Frailty model: No covariates \n")
			cat("    -------------------------- \n")
			cat("\n")
		}
		
		cat("\n")
		
		cat("    Convergence criteria: \n")
		cat("    parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
		
		cat("\n")
		cat("    n=", x$n)
		if (length(x$na.action)){ 
			cat("  (", length(x$na.action), " observation deleted due to missing) \n")
		}else{
			cat("\n")
		}
		cat("    n events=", x$n.events, " n groups=", x$groups)
		cat( "\n")
		cat("    number of iterations: ", x$n.iter)	
	}
	cat("\n")
	invisible()
}
