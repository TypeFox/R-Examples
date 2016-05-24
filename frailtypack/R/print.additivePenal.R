
"print.additivePenal" <- function (x, digits = max(options()$digits - 4, 6), ...) 
{
	if (!is.null(cl <- x$call)){
		cat("Call:\n")
		dput(cl)
		if (x$type == "counting" & x$AG == FALSE){
			cat("\n      left truncated structure used")
		}
		if (x$AG == TRUE){
			cat("\n      Calendar timescale")
		}
		cat("\n")
	}
	if (!is.null(x$fail)){
		cat(" additivePenal failed.", x$fail, "\n")
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
			cat("\n")
		}  
		if (x$n.knots.temp > 20){
			cat("\n")
			cat("  The maximum number of knots is 20","\n")
		} 
	}else{
		if ((x$typeof == 1) & (x$indic.nb.int == 1)) cat("  The maximum number of time intervals is 20","\n")
	}
#AD
 	if (x$istop == 1){
		if (!is.null(coef)){ 
# JRG sep '09
#        if (nvar>1)
#         {    
#          seH <- sqrt(diag(x$varH))
#          seHIH <- sqrt(diag(x$varHIH))
#         }
#        if (nvar==1)
#         {    
			seH <- sqrt(x$varH)
			seHIH <- sqrt(x$varHIH)
#         }
			if (x$typeof == 0){
				tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - 
				pchisq((coef/seH)^2, 1), digits - 1))
				if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,x$p.global_chisq)
			}else{
				tmp <- cbind(coef, exp(coef), seH, coef/seH, signif(1 - 
				pchisq((coef/seH)^2, 1), digits - 1))
				if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,x$p.global_chisq)
			}
	
			cat("\n")
			cat("  Additive gaussian frailty model parameter estimates ","\n")	
			if (x$typeof == 0){
				cat("  using a Penalized Likelihood on the hazard function","\n")
			}else{	
				cat("  using a Parametrical approach for the hazard function","\n")
			}
	
			if (x$n.strat>1) cat("  (Stratification structure used)", "\n")
			if (x$typeof == 0){
				if(x$global_chisq.test==1){
					dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
					
				}
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "SE coef (HIH)", "z", "p"))

			}else{
				if(x$global_chisq.test==1){
					dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
					
				}
				dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","SE coef (H)", "z", "p"))

			}
			cat("\n")
			prmatrix(tmp)
			if(x$global_chisq.test==1){
				cat("\n")
				prmatrix(tmpwald)
			}
			cat("\n")


			if (x$correlation){
					cat("    Covariance (between the two frailty terms,  \n")
					cat("        the intercept and the slope):", x$cov, "(SE:",x$varcov^0.5, ")", "\n")
					cat("\n")
			}
		} 
		
#AD:		
		if (x$noVar == 1){
			cat("\n")
			cat("    Additive gaussian frailty model: No covariates \n")
			cat("    ------------------------------- \n")
			cat("\n")
		}
#AD:
		if (x$rho!=-1)cat("    Corresponding correlation between the two frailty terms :", x$rho, "\n")
		if (x$typeof == 0){
			cat("    Variance for random intercept:", x$sigma2, "(SE (H):", x$varSigma2[1]^.5, ")", "\n") # "(SE (HIH):", x$varSigma2[2]^.5, ")", "\n")
	
			cat("    Variance for random slope:", x$tau2, "(SE (H):", x$varTau2[1]^.5, ")", "\n") #"(SE (HIH):", x$varTau2[2]^.5, ")", "\n")
		}else{
			cat("    Variance for random intercept:", x$sigma2, "(SE (H):", x$varSigma2[1]^.5, ")", "\n")
	
			cat("    Variance for random slope:", x$tau2, "(SE (H):", x$varTau2[1]^.5, ")", "\n")
		}
	
	
	
		cat(" \n")
		if (x$typeof == 0){
			cat(paste("    penalized marginal log-likelihood =", round(x$logLikPenal,2)))
#AD:
			cat("\n")
			cat("    Convergence criteria: \n")
			cat("    parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
			cat("\n")
			cat("    LCV = the approximate likelihood cross-validation criterion\n")
			cat("    in the semi parametrical case     =",x$LCV,"\n")
#AD:
		}else{
			cat(paste("    marginal log-likelihood =", round(x$logLik,2)))
			cat("\n")
			cat("      Convergence criteria: \n")
			cat("      parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
			cat("\n")
#			cat("      LCV = the approximate likelihood cross-validation criterion\n")
#			cat("            in the parametrical case     =",x$LCV,"\n")	
			cat("      AIC = Aikaike information Criterion     =",x$AIC,"\n")
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
		cat("\n")
		cat("    n=", x$n)
		if (length(x$na.action)){ 
			cat("  (", length(x$na.action), " observation deleted due to missing) \n")
		}else{ 
			cat("\n")
		}
		cat("    n events=", x$n.event, " n groups=", x$groups)
		cat( "\n")
		cat("    number of iterations: ", x$n.iter,"\n")
		if ((x$typeof == 1) & (x$indic.nb.int == 1)){
			cat("      Exact number of time intervals used: 20","\n")
		 }else{
		 	if (x$typeof == 1) cat("      Exact number of time intervals used: ",x$nbintervR,"\n")
		 }
		cat("\n")
		if (x$typeof == 0){
			cat("    Exact number of knots used: ", x$n.knots, "\n")
	
			if (!x$cross.Val){
				cat("    Value of the smoothing parameter: ", x$kappa, sep=" ")
			}
	
			if (x$cross.Val){
				if (is.null(x$theta)){
					cat("    Smoothing parameter estimated by Cross validation: ", x$kappa, sep=" ")
				}else{
					cat("    Best smoothing parameter estimated by")
					cat("\n")
					cat("       an approximated Cross validation: ", x$kappa, sep=" ")
				} 
			}
		cat(", DoF: ", formatC(-x$DoF, format="f",digits=2))
		}
	}else{
		if (!is.null(coef)){
			cat("\n")
			cat("  Additive gaussian frailty model parameter estimates ","\n")	
			if (x$typeof == 0){
				cat("  using a Penalized Likelihood on the hazard function","\n")
			}else{	
				cat("  using a Parametrical approach for the hazard function","\n")	
			}
	
			if (x$n.strat>1) cat("  (Stratification structure used)", "\n")    	
		}
		if (x$noVar == 1){
			cat("\n")
			cat("    Additive gaussian frailty model: No covariates \n")
			cat("    ------------------------------- \n")
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
		cat("    n events=", x$n.event, " n groups=", x$groups)
		cat( "\n")
		cat("    number of iterations: ", x$n.iter)
		cat("\n")
	}
	
	
	cat("\n")
	invisible()
}
