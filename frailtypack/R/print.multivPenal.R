
"print.multivPenal" <- function (x, digits = max(options()$digits - 4, 6), ...) 
{
	if (!is.null(cl <- x$call)){
		cat("Call:\n")
		dput(cl)
#  		if (x$type == "counting"){
#  			cat("\n      left truncated structure used")
#  		}
		if (x$AG == TRUE){
			cat("\n      Calendar timescale")
		}
		cat("\n")
	}
	if (!is.null(x$fail)) {
		cat(" multivPenal failed.", x$fail, "\n")
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
		if (all(x$n.knots.temp < 4)){
			cat("\n")
			cat("  The minimum number of knots is 4","\n")	
			cat("\n")
		} 
		if (all(x$n.knots.temp > 20)){
			cat("\n")
			cat("  The maximum number of knots is 20","\n")	
		}  
	}else{
		if ((x$typeof == 1) & (x$indic.nb.int1 == 1)) cat("  The maximum number of time intervals nb.int[1] is 20","\n")
		if ((x$typeof == 1) & (x$indic.nb.int2 == 1)) cat("  The maximum number of time intervals nb.int[2] is 20","\n")
		if ((x$typeof == 1) & (x$indic.nb.int3 == 1)) cat("  The maximum number of time intervals nb.int[3] is 20","\n")
	}
#AD
	if (x$istop == 1){
		if (!is.null(coef)){ 
			seH <- sqrt(diag(x$varH))[-c(1,2)]
			seHIH <- sqrt(diag(x$varHIH))[-c(1,2)]
			if (x$typeof == 0){
				tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - 
				pchisq((coef/seH)^2, 1), digits - 1))
				if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,x$p.global_chisq)
				if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d,x$dof_chisq_d,x$p.global_chisq_d)
				if(x$global_chisq.test2==1) tmpwald2 <- cbind(x$global_chisq2,x$dof_chisq2,x$p.global_chisq2)
			}else{
				tmp <- cbind(coef, exp(coef), seH, coef/seH, signif(1 - 
				pchisq((coef/seH)^2, 1), digits - 1))
				if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,x$p.global_chisq)
				if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d,x$dof_chisq_d,x$p.global_chisq_d)
				if(x$global_chisq.test2==1) tmpwald2 <- cbind(x$global_chisq2,x$dof_chisq2,x$p.global_chisq2)
			}
			cat("\n")
			cat("  Multivariate joint gaussian frailty model for two survival outcomes and a terminal event","\n")
			if (x$typeof == 0){
				cat("  using a Penalized Likelihood on the hazard function","\n")
			}else{
				cat("  using a Parametrical approach for the hazard function","\n")
			}
			if (x$typeof == 0){
				if(x$global_chisq.test==1){
					dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
					
				}
				if(x$global_chisq.test_d==1){
					dimnames(tmpwalddc) <- list(x$names.factordc,c("chisq", "df", "global p"))
					
				}
				if(x$global_chisq.test2==1){
					dimnames(tmpwald2) <- list(x$names.factor2,c("chisq", "df", "global p"))
					
				}

				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
				"SE coef (H)", "SE coef (HIH)", "z", "p"))
			}else{
				if(x$global_chisq.test==1){
					dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
					
				}
				if(x$global_chisq.test_d==1){
					dimnames(tmpwalddc) <- list(x$names.factordc,c("chisq", "df", "global p"))
					
				}
				if(x$global_chisq.test2==1){
					dimnames(tmpwald2) <- list(x$names.factor2,c("chisq", "df", "global p"))
					
				}
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", 
				"SE coef (H)", "z", "p"))
			}
	
			if (x$noVar[1] == 0){
				cat("\n")
				cat("Recurrences 1:\n")
				cat("------------ \n")
				prmatrix(tmp[1:x$nvar[1], ,drop=FALSE])
				if(x$global_chisq.test==1){
					cat("\n")
					prmatrix(tmpwald)
				}
			}
			cat("\n")

			if (x$noVar[3] == 0){	
				cat("Recurrences 2:\n")
				cat("------------- \n")
				ntmp <- sum(x$nvar[1:2])
				prmatrix(tmp[-c(1:ntmp), ,drop=FALSE])
				if(x$global_chisq.test2==1){
					cat("\n")
					prmatrix(tmpwald2)
				}
			cat("\n")

			if (x$noVar[2] == 0){	
				cat("Terminal event:\n")
				cat("---------------- \n")
				ntmp <- sum(x$nvar[1:2])
				ntmp2 <- sum(x$nvar)
				prmatrix(tmp[-c(1:x$nvar[1],(ntmp+1):ntmp2), ,drop=FALSE])
				if(x$global_chisq.test_d==1){
					cat("\n")
					prmatrix(tmpwalddc)
				}
				cat("\n")
			}

			}
		}

#theta
		theta1 <- x$theta1
		temp <- diag(x$varH)[1]
		seH.theta1 <- sqrt(((2 * (theta1^0.5))^2) * temp)
		temp <- diag(x$varHIH)[1]
		seHIH.theta1 <- sqrt(((2 * (theta1^0.5))^2) * temp)

#eta
		theta2 <- x$theta2
		temp <- diag(x$varH)[2]
		seH.theta2 <- sqrt(((2 * (theta2^0.5))^2) * temp)
		temp <- diag(x$varHIH)[2]
		seHIH.theta2 <- sqrt(((2 * (theta2^0.5))^2) * temp)

#AD:
		if (x$noVar[1] == 1){
			cat("\n")
			cat("    Recurrences: No covariates \n")
			cat("    ----------- \n")
		}
		
		if (x$noVar[2] == 1){
			cat("\n")
			cat("    Terminal event: No covariates \n")
			cat("    -------------- \n")
			cat("\n")
		}
		if (x$noVar[3] == 1){
			cat("\n")
			cat("    Recurrences2: No covariates \n")
			cat("    ----------- \n")
		}
#AD:  
		cat(" Parameters associated with Frailties: \n")
# 		if (x$typeof == 0){
# 			cat("   theta1 :", theta1, "(SE (H):", seH.theta1, ")", "(SE (HIH):", seHIH.theta1, ")", "\n")
# 			cat("   theta2 :", theta2, "(SE (H):", seH.theta2, ")", "(SE (HIH):", seHIH.theta2, ")", "\n")
# 			cat("   alpha1 :", x$alpha1, "(SE (H):", sqrt(diag(x$varH))[3], ")", "(SE (HIH):", sqrt(diag(x$varHIH))[3],")", "\n")
# 			cat("   alpha2 :", x$alpha2, "(SE (H):", sqrt(diag(x$varH))[4], ")", "(SE (HIH):", sqrt(diag(x$varHIH))[4],")", "\n")
# 			cat("   rho :", x$rho, "(SE (H):", sqrt(diag(x$varH))[5], ")", "(SE (HIH):", sqrt(diag(x$varHIH))[5],")", "\n")
# 			cat(" \n")
# 		}else{
# 			cat("   theta1 :", theta1, "(SE (H):", seH.theta1, ")", "\n")
# 			cat("   theta2 :", theta2, "(SE (H):", seH.theta2, ")", "\n")
# 			cat("   alpha1 :", x$alpha1, "(SE (H):", sqrt(diag(x$varH))[3], ")", "\n")
# 			cat("   alpha2 :", x$alpha2, "(SE (H):", sqrt(diag(x$varH))[4], ")", "\n")
# 			cat("   rho :", x$rho, "(SE (H):", sqrt(diag(x$varH))[5], ")", "\n")
# 			cat(" \n")
# 		}
		cat("   theta1 :", theta1, "(SE (H):", seH.theta1, ")", "p =", signif(1 - pnorm(theta1/seH.theta1), digits - 1), "\n")
		cat("   theta2 :", theta2, "(SE (H):", seH.theta2, ")", "p =", signif(1 - pnorm(theta2/seH.theta2), digits - 1), "\n")
		cat("   alpha1 :", x$alpha1, "(SE (H):", sqrt(diag(x$varH))[3], ")", "p =", signif(1 - pchisq((x$alpha1/sqrt(diag(x$varH))[3])^2,1), digits - 1), "\n")
		cat("   alpha2 :", x$alpha2, "(SE (H):", sqrt(diag(x$varH))[4], ")", "p =", signif(1 - pchisq((x$alpha2/sqrt(diag(x$varH))[4])^2,1), digits - 1), "\n")
		cat("   rho :", x$rho, "(SE (H):", sqrt(diag(x$varH))[5], ")", "\n")
		cat(" \n")

		if (x$typeof == 0){
			cat(paste("   penalized marginal log-likelihood =", round(x$logLikPenal,2)))
			cat("\n")
			cat("   LCV = the approximate likelihood cross-validation criterion\n")
			cat("         in the semi parametric case     =",x$LCV,"\n")
		}else{
			cat(paste("   marginal log-likelihood =", round(x$logLik,2)))
			cat("\n")
#			cat("   LCV = the approximate likelihood cross-validation criterion\n")
#			cat("         in the parametric case     =",x$LCV,"\n")	
			cat("   AIC = Aikaike information Criterion     =",x$AIC,"\n")
			cat("\n")
			cat("The expression of the Aikaike Criterion is:","\n")
			cat("        'AIC = (1/n)[np - l(.)]'","\n")
			if (x$typeof == 2){
				cat("\n")
				cat("      Scale for the weibull hazard function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),round(x$scale.weib[3],2),"\n")	
				cat("      Shape for the weibull hazard function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),round(x$shape.weib[3],2),"\n")
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
		cat("   n=", x$n)
		if (length(x$na.action)){
			cat("      (", length(x$na.action), " observation deleted due to missing) \n")
		}else{ 
			cat("\n")
		}
		cat("   n recurrent events of type 1=", x$n.events, " n subjects=", x$groups)
		cat("\n")
		cat("   n recurrent events of type 2=", x$n.events2)
		cat("\n")
		cat("   n terminal events=", x$n.deaths)
		cat("\n")
		cat("   number of iterations: ", x$n.iter,"\n")
		
		if ((x$typeof == 1) & (x$indic.nb.int1 == 1)){
			cat("   Exact number of intervals recurrent of type 1 used for hazard: 20","\n")
		 }else{
		 	if (x$typeof == 1) cat("   Exact number of intervals for recurrent of type 1 used for hazard: ",x$nbintervR,"\n")
		 }

		if ((x$typeof == 1) & (x$indic.nb.int3 == 1)){
			cat("   Exact number of intervals for recurrent of type 2 used for hazard: 20","\n")
		 }else{
		 	if (x$typeof == 1) cat("   Exact number of intervals for recurrent of type 2 used for hazard: ",x$nbintervR2,"\n")
		 }

		if ((x$typeof == 1) & (x$indic.nb.int2 == 1)){
			cat("   Exact number of intervals for terminal event used for hazard: 20","\n")
		 }else{
		 	if (x$typeof == 1) cat("   Exact number of intervals for terminal event used for hazard: ",x$nbintervDC,"\n")
		 }
		

		if (x$typeof == 0){ 
			cat("\n")
			cat("   Exact number of knots used: ", x$n.knots[1]," ",x$n.knots[2]," ",x$n.knots[3], "\n")
			cat("   Value of the smoothing parameters: ", x$kappa, sep=" ")
			cat("\n")
		}
	}else{
		if (!is.null(coef)){ 
			cat("\n")
			cat("  Multivariate joint gaussian frailty model for two survival outcomes and a terminal event","\n")
			if (x$typeof == 0){
				cat("  using a Penalized Likelihood on the hazard function","\n")
			}else{
				cat("  using a Parametrical approach for the hazard function","\n")
			}	
			
			if (x$noVar[1] == 1){
				cat("\n")
				cat("    Recurrences: No covariates \n")
				cat("    ----------- \n")
			}
			
			if (x$noVar[3] == 1){
				cat("\n")
				cat("    Recurrences 2: No covariates \n")
				cat("    -------------- \n")
				cat("\n")
			}

			if (x$noVar[2] == 1){
				cat("\n")
				cat("    Terminal event: No covariates \n")
				cat("    -------------- \n")
				cat("\n")
			}

			cat("\n")
			cat("   n=", x$n)
			if (length(x$na.action)){
				cat("      (", length(x$na.action), " observation deleted due to missing) \n")
			}else{ 
				cat("\n")
			}
			cat("   n recurrent events of type 1=", x$n.events, " n subjects=", x$groups)
			cat("\n")
			cat("   n recurrent events of type 2 =", x$n.events2)
			cat("\n")
			cat("   n terminal events=", x$n.deaths)
			cat("\n")
		}
	}
    invisible()
}
