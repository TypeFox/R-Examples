
"print.jointPenal" <- function (x, digits = max(options()$digits - 4, 6), ...) 
{
  
  if (x$istop == 1){
    # plot des coefficient dependant du temps
    if (any(x$nvartimedep != 0)) par(mfrow=c(1,2))#sum(as.logical(x$nvartimedep)),max(x$nvartimedep)))
    if ((x$nvartimedep[1] != 0) & (x$istop == 1)){
      for (i in 0:(x$nvartimedep[1]-1)){
        matplot(x$BetaTpsMat[,1],x$BetaTpsMat[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=paste("Recurrent : ",x$Names.vardep[i+1]),ylim=c(min(x$BetaTpsMat[,-1]),max(x$BetaTpsMat[,-1])))
      }
    }
    if ((x$nvartimedep[2] != 0) & (x$istop == 1)){
      for (i in 0:(x$nvartimedep[2]-1)){
        matplot(x$BetaTpsMatDc[,1],x$BetaTpsMatDc[,(2:4)+4*i],col="blue",type="l",lty=c(1,2,2),xlab="t",ylab="beta(t)",main=paste("Death : ",x$Names.vardepdc[i+1]),ylim=c(min(x$BetaTpsMatDc[,-1]),max(x$BetaTpsMatDc[,-1])))
      }
    }
  }
  
  if (!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl)
    #if (x$type == "counting" & x$AG == FALSE){ # pour l'instant joint n'accepte pas la vraie troncature a gauche
    #	cat("\n      left truncated structure used")
    #}
    if (x$AG == TRUE){
      cat("\n      Calendar timescale")
    }
    if (x$intcens == TRUE){
      cat("\n      interval censored data used")
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
      cat("\n")
    } 
    if (x$n.knots.temp > 20){
      cat("\n")
      cat("  The maximum number of knots is 20","\n")
    }
  }else{
    if ((x$typeof == 1) & (x$indic.nb.intR == 1)) cat("  The maximum number of time intervals is 20","\n")
    if ((x$typeof == 1) & (x$indic.nb.intD == 1)) cat("  The maximum number of time intervals is 20","\n")
  }
  #AD
  
  if (x$logNormal == 0) frail <- x$theta
  else frail <- x$sigma2
  indic_alpha <- x$indic_alpha
  
  if (x$istop == 1){
    if (!is.null(coef)){
      if (indic_alpha == 1 || x$joint.clust==2) {
        seH <- sqrt(diag(x$varH))[-c(1,2)]
        seHIH <- sqrt(diag(x$varHIH))[-c(1,2)]
      }else{
        seH <- sqrt(diag(x$varH))[-1]
        seHIH <- sqrt(diag(x$varHIH))[-1]
      }
      if (x$typeof == 0){
        tmp <- cbind(coef, exp(coef), seH, seHIH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits - 1))
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,x$p.global_chisq)
        if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d,x$dof_chisq_d,x$p.global_chisq_d)
      }else{
        tmp <- cbind(coef, exp(coef), seH, coef/seH, signif(1 - pchisq((coef/seH)^2, 1), digits - 1))
        if(x$global_chisq.test==1) tmpwald <- cbind(x$global_chisq,x$dof_chisq,x$p.global_chisq)
        if(x$global_chisq.test_d==1) tmpwalddc <- cbind(x$global_chisq_d,x$dof_chisq_d,x$p.global_chisq_d)
      }
      cat("\n")
      if (x$joint.clust == 0) cat("  For clustered data","\n")
      if (x$joint.clust == 0){
        if (x$logNormal == 0){
          cat("  Joint gamma frailty model for a survival and a terminal event processes","\n")
        }else{
          cat("  Joint Log-Normal frailty model for a survival and a terminal event processes","\n")
        }
      }else{
        if ((x$logNormal == 0)&(x$joint.clust==1)){
          cat("  Joint gamma frailty model for recurrent and a terminal event processes","\n")
        }
        else if ((x$logNormal == 0)&(x$joint.clust==2)){
          cat("  General Joint gamma frailty model for recurrent and a terminal event processes","\n")
        }
        else{
          cat("  Joint Log-Normal frailty model for recurrent and a terminal event processes","\n")
        }
      }
      if (x$typeof == 0){
        cat("  using a Penalized Likelihood on the hazard function","\n")
      }else{
        cat("  using a Parametrical approach for the hazard function","\n")
      }
      if (any(x$nvartimedep != 0)) cat("  and some time-dependant covariates","\n")
      if (x$n.strat>1) cat("  (Stratification structure used for recurrences) :",x$n.strat,"strata \n")
      if (x$typeof == 0){
        if(x$global_chisq.test==1){
          dimnames(tmpwald) <- list(x$names.factor,c("chisq", "df", "global p"))
          
        }
        if(x$global_chisq.test_d==1){
          dimnames(tmpwalddc) <- list(x$names.factordc,c("chisq", "df", "global p"))
          
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
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                             "SE coef (H)", "z", "p"))
      }
      cat("\n")
      
      if (x$nvarnotdep[1] == 0){
        if (x$joint.clust == 0) cat("Survival event:\n")
        if ((x$joint.clust == 1) | (x$joint.clust == 2)) cat("Recurrences:\n")
        cat("------------- \n")
        cat("No constant coefficients, only time-varying effects of the covariates \n")
      }else{
        if (x$noVar1 == 0){
          if (x$joint.clust == 0) cat("Survival event:\n")
          if (x$joint.clust >= 1) cat("Recurrences:\n")
          cat("------------- \n")
          prmatrix(tmp[1:x$nvarnotdep[1], ,drop=FALSE])
          if(x$global_chisq.test==1){
            cat("\n")
            prmatrix(tmpwald)
          }
        }
      }
      cat("\n")
      
      if (x$nvarnotdep[2] == 0){
        cat("Terminal event:\n")
        cat("---------------- \n")
        cat("No constant coefficients, only time-varying effects of the covariates \n")
      }else{
        if (x$noVar2 == 0){
          cat("Terminal event:\n")
          cat("---------------- \n")
          prmatrix(tmp[-c(1:x$nvarnotdep[1]), ,drop=FALSE])
          if(x$global_chisq.test_d==1){
            cat("\n")
            prmatrix(tmpwalddc)
          }
        }
      }
      cat("\n")
    }
    #theta <- x$theta
    temp <- diag(x$varH)[1]
    seH.frail <- sqrt(((2 * (frail^0.5))^2) * temp) # delta methode
    temp <- diag(x$varHIH)[1]
    seHIH.frail <- sqrt(((2 * (frail^0.5))^2) * temp) # delta methode
    #AD:
    if (x$noVar1 == 1){
      cat("\n")
      if (x$joint.clust == 0) cat("    Survival event: No covariates \n")
      if (x$joint.clust >= 1) cat("    Recurrences: No covariates \n")
      cat("    ----------- \n")
    }
    
    if (x$noVar2 == 1){
      cat("\n")
      cat("    Terminal event: No covariates \n")
      cat("    -------------- \n")
      cat("\n")
    }
    #AD:  
    cat(" Frailty parameters: \n")
    # 		if (x$typeof == 0){
    # 			if (x$logNormal == 0){
    # 				cat("   theta (variance of Frailties, w):", frail, "(SE (H):",seH.frail, ")", "(SE (HIH):", seHIH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (w^alpha for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")","(SE (HIH):",sqrt(diag(x$varHIH))[2],")","\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}else{
    # 				cat("   sigma square (variance of Frailties, eta):", frail, "(SE (H):",seH.frail, ")", "(SE (HIH):", seHIH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (exp(alpha.eta) for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")","(SE (HIH):",sqrt(diag(x$varHIH))[2],")","\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}
    # 			cat(" \n")
    # 		}else{
    # 			if (x$logNormal == 0){
    # 				cat("   theta (variance of Frailties, Z):", frail, "(SE (H):",seH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (w^alpha for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")", "\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}else{
    # 				cat("   sigma square (variance of Frailties, Z):", frail, "(SE (H):",seH.frail, ")", "\n")
    # 				if (indic_alpha == 1) cat("   alpha (exp(alpha.eta) for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")", "\n")
    # 				else cat("   alpha is fixed (=1) \n")
    # 			}
    # 			cat(" \n")
    # 		}
    if (x$logNormal == 0){
      cat("   theta (variance of Frailties, w):", frail, "(SE (H):",seH.frail, ")", "p =", signif(1 - pnorm(frail/seH.frail), digits - 1), "\n")
      if (indic_alpha == 1 & x$joint.clust<=1) cat("   alpha (w^alpha for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")", "p =", signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[2])^2,1), digits - 1), "\n")
      else if (x$joint.clust ==2) cat("   eta   (variance of Frailties, v):", x$eta, "(SE (H):",sqrt(((2 * (x$eta^0.5))^2) * diag(x$varH)[2]), ")", "p =", signif(1 - pnorm (x$eta/sqrt(((2 * (x$eta^0.5))^2) * diag(x$varH)[2]),1), digits - 1), "\n")
      else cat("   alpha is fixed (=1) \n")
    }else{
      cat("   sigma square (variance of Frailties, eta):", frail, "(SE (H):",seH.frail, ")", "p =", signif(1 - pnorm(frail/seH.frail), digits - 1), "\n")
      if (indic_alpha == 1) cat("   alpha (exp(alpha.eta) for terminal event):", x$alpha, "(SE (H):",sqrt(diag(x$varH))[2], ")", "p =", signif(1 - pchisq((x$alpha/sqrt(diag(x$varH))[2])^2,1), digits - 1), "\n")
      else cat("   alpha is fixed (=1) \n")
    }
    cat(" \n")
    
    if (x$typeof == 0){
      cat(paste("   penalized marginal log-likelihood =", round(x$logLikPenal,2)))
      cat("\n")
      cat("   Convergence criteria: \n")
      cat("   parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
      cat("\n")
      cat("   LCV = the approximate likelihood cross-validation criterion\n")
      cat("         in the semi parametric case     =",x$LCV,"\n")
    }else{
      cat(paste("   marginal log-likelihood =", round(x$logLik,2)))
      cat("\n")
      cat("   Convergence criteria: \n")
      cat("   parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
      cat("\n")
      #			cat("   LCV = the approximate likelihood cross-validation criterion\n")
      #			cat("         in the parametric case     =",x$LCV,"\n")
      cat("   AIC = Aikaike information Criterion     =",x$AIC,"\n")
      cat("\n")
      cat("The expression of the Aikaike Criterion is:","\n")
      cat("        'AIC = (1/n)[np - l(.)]'","\n")
      if (x$typeof == 2){
        cat("\n")
        cat("      Scale for the weibull hazard function is :",round(x$scale.weib[1],2),round(x$scale.weib[2],2),"\n")	
        cat("      Shape for the weibull hazard function is :",round(x$shape.weib[1],2),round(x$shape.weib[2],2),"\n")
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
    if (x$joint.clust == 0){
      cat("   n observations=", x$n, " n subjects=", x$ind, " n groups=", x$groups)
    }else{
      cat("   n observations=", x$n, " n subjects=", x$groups)
    }
    if (length(x$na.action)){
      cat("      (", length(x$na.action), " observation deleted due to missing) \n")
    }else{ 
      cat("\n")
    }
    if (x$joint.clust == 0){
      cat("   n events=", x$n.events)
    }else{
      cat("   n recurrent events=", x$n.events)
    }
    cat("\n")
    cat("   n terminal events=", x$n.deaths)
    cat("\n")
    cat("   n censored events=" ,x$n.censored)
    cat("\n")
    cat("   number of iterations: ", x$n.iter,"\n")
    
    if ((x$typeof == 1) & (x$indic.nb.intR == 1)){
      cat("   Exact number of time intervals used: 20","\n")
    }else{
      if (x$typeof == 1) cat("   Exact number of time intervals used: ",x$nbintervR,"\n")
    }
    if ((x$typeof == 1) & (x$indic.nb.intD == 1)){
      cat("   Exact number of time intervals used: 20","\n")
    }else{
      if (x$typeof == 1) cat("   Exact number of time intervals used: ",x$nbintervDC,"\n")
    }
    
    if (x$typeof == 0){ 
      cat("\n")
      cat("   Exact number of knots used: ", x$n.knots, "\n")
      cat("   Value of the smoothing parameters: ", x$kappa, sep=" ")
      cat("\n")
    }
  }else{
    if (!is.null(coef)){ 
      cat("\n")
      if (x$joint.clust == 0) cat("  For clustered data","\n")
      if (x$joint.clust == 0){
        if (x$logNormal == 0){
          cat("  Joint gamma frailty model for a survival and a terminal event processes","\n")
        }else{
          cat("  Joint Log-Normal frailty model for a survival and a terminal event processes","\n")
        }
      }else{
        if ((x$logNormal == 0)&(x$joint.clust==1)){
          cat("  Joint gamma frailty model for recurrent and a terminal event processes","\n")
        }
        else if ((x$logNormal == 0)&(x$joint.clust==2)){
          cat("  General Joint gamma frailty model for recurrent and a terminal event processes","\n")
        }
        else{
          cat("  Joint Log-Normal frailty model for recurrent and a terminal event processes","\n")
        }
      }
      if (x$typeof == 0){
        cat("  using a Penalized Likelihood on the hazard function","\n")
      }else{
        cat("  using a Parametrical approach for the hazard function","\n")
      }
      if (any(x$nvartimedep != 0)) cat("  and some time-dependant covariates","\n")
      if (x$noVar1 == 1){
        cat("\n")
        if (x$joint.clust == 0) cat("    Survival event: No covariates \n")
        if (x$joint.clust >= 1) cat("    Recurrences: No covariates \n")
        cat("    ----------- \n")
      }
      
      if (x$noVar2 == 1){
        cat("\n")
        cat("    Terminal event: No covariates \n")
        cat("    -------------- \n")
        cat("\n")
      }
      
      cat("\n")
      
      cat("   Convergence criteria: \n")
      cat("   parameters =",signif(x$EPS[1],3),"likelihood =",signif(x$EPS[2],3),"gradient =",signif(x$EPS[3],3),"\n")
      
      cat("\n")
      cat("   n=", x$n)
      if (length(x$na.action)){
        cat("      (", length(x$na.action), " observation deleted due to missing) \n")
      }else{
        cat("\n")
      }
      if (x$joint.clust == 0){
        cat("   n events=", x$n.events)
      }else{
        cat("   n recurrent events=", x$n.events)
      }
      cat("\n")
      cat("   n terminal events=", x$n.deaths)
      cat("\n")
    }
  }
  invisible()
}
