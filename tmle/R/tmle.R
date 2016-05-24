# July 3, 2012
# With Noemi - add starting values for fitting epsilon start = c(0,0)
# Targeted Maximum Likelihood Estimation
# for non-parametric estimation of the marginal effect of a binary point
# treatment, adjusting for treatment (g) and missingness (g.Delta) mechanisms
# Parameters include: additive treatment effect: E_W[E(Y|A=1,W) - E(Y|A=0,W)]
# and, for binary outcomes, relative risk(RR) and odds ratio(RR)
#   mu1 = E_W[E(Y|A=1,W)], mu0 = E_W[E(Y|A=0,W)]
#   RR = mu1/mu0
#   OR =  mu1/(1-mu1)/(mu0 / (1-mu0))
# Controlled direct effect estimation available for optional binary Z intermediate variable,
# (P(DeltaY|Z,A,W, Delta=1)P(Delta|Z,A,W)P(Z|A,W)P(A|W)P(W))
# EY1 parameter estimated when there is missinginess and no treatment assignment
# author: Susan Gruber, sgruber@berkeley.edu
# date:   September 25, 2010
# revised: December 20, 2010
#          March 22, 2012 
# Models or estimated values for Q, g, g.Z, g.Delta can be user-supplied or 
# estimated using super learner.  Cross-validated inital Q can be obtained
# by specifying cvQinit=TRUE

tmleNews <- function(...){
	RShowDoc("NEWS", package="tmle",...)
}

#-------------summary.tmle------------------
# create tmle summary object
# return values of psi, variance, epsilon, 
# and info on estimation of Q and g factors, if available
#---------------------------------------	 
summary.tmle <- function(object,...) {
  if(identical(class(object), "tmle")){
		Qmodel <- Qcoef <- gmodel <- gcoef <- g.Deltamodel <- g.Deltacoef <- NULL
		g.Zmodel <- g.Zcoef <- NULL
		Qterms <- gterms <- g.Deltaterms <- g.Zterms <- ""
		if(!is.null(object$Qinit)){
		  if (!is.null(object$Qinit$coef)) {
			Qcoef <- object$Qinit$coef	
			if(class(Qcoef) == "matrix"){
				Qterms <- colnames(Qcoef)
			} else {
				Qterms <- names(Qcoef)
			}
			Qmodel <- paste("Y ~ 1")
			if(length(Qterms) > 1) {
				Qmodel <- paste("Y ~ ", paste(Qterms, collapse =" + "))
			} 
		  } else {
		  	Qmodel <- object$Qinit$type
		  }
		}
		if(!is.null(object$g)){
		if (!is.null(object$g$coef)) {
			gbd <- object$g$bound
			gcoef <- object$g$coef
			if(class(gcoef) == "matrix"){
				gterms <- colnames(gcoef)
			} else {
				gterms <- names(gcoef)
			}			
			gmodel <- paste("A ~ 1")
			if(length(gterms) > 1) {
				gmodel <- paste("A ~ ", paste(gterms, collapse =" + "))
			} 
		}}
		if(!is.null(object$g.Z)){
		if (!is.null(object$g.Z$coef)) {
			g.Zcoef <- object$g.Z$coef
			if(class(g.Zcoef) == "matrix"){
				g.Zterms <- colnames(g.Zcoef)
			} else {
				g.Zterms <- names(g.Zcoef)
			}			
			g.Zmodel <- paste("Z ~ 1")
			if(length(g.Zterms) > 1) {
				g.Zmodel <- paste("Z ~", paste(g.Zterms, collapse =" + "))
			} 
		}}
		if(!is.null(object$g.Delta)){
		if (!is.null(object$g.Delta$coef)) {
			g.Deltacoef <- object$g.Delta$coef
			if(class(g.Deltacoef) == "matrix"){
				g.Deltaterms <- colnames(g.Deltacoef)
			} else {
				g.Deltaterms <- names(g.Deltacoef)
			}			
			g.Deltamodel <- paste("Delta ~ 1")
			if(length(g.Deltaterms) > 1) {
				g.Deltamodel <- paste("Delta ~", paste(g.Deltaterms, collapse =" + "))
			} 
		}}
		summary.tmle <- list(estimates=object$estimates,
						Qmodel=Qmodel, Qterms=Qterms, Qcoef=Qcoef, Qtype=object$Qinit$type,
						gbd=gbd, gmodel=gmodel, gterms=gterms, gcoef=gcoef,gtype=object$g$type, 
						g.Zmodel=g.Zmodel, g.Zterms=g.Zterms, g.Zcoef=g.Zcoef, g.Ztype=object$g.Z$type,
						g.Deltamodel=g.Deltamodel, g.Deltaterms=g.Deltaterms, g.Deltacoef=g.Deltacoef, g.Deltatype=object$g.Delta$type )
		class(summary.tmle) <- "summary.tmle"
	} else {
		stop("object must have class 'tmle'")
		summary.tmle <- NULL
	}
	return(summary.tmle)
}

#-------------print.summary.tmle------------------
# print tmle summary object
#-------------------------------------------------
print.summary.tmle <- function(x,...) {
  if(identical(class(x), "summary.tmle")){
  	cat(" Initial estimation of Q\n")
  	cat("\t Procedure:", x$Qtype)
  	if(!(is.na(x$Qcoef[1]))){
  		cat("\n\t Model:\n\t\t", x$Qmodel)
  		cat("\n\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$Qterms)
  		extra <- ifelse(x$Qcoef >= 0, "  ", " ")
  		for(i in 1:length(x$Qcoef)){
  			cat("\t", terms[i], extra[i], x$Qcoef[i], "\n")
  		}
  	}
  	cat("\n Estimation of g (treatment mechanism)\n")
  	cat("\t Procedure:", x$gtype,"\n")
  	if(!(is.na(x$gcoef[1]))){	
 		cat("\t Model:\n\t\t", x$gmodel, "\n")
  		cat("\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$gterms)
  		extra <- ifelse(x$gcoef >= 0, "  ", " ")
  		for(i in 1:length(x$gcoef)){
  			cat("\t", terms[i], extra[i], x$gcoef[i], "\n")
  	}}
  	cat("\n Estimation of g.Z (intermediate variable assignment mechanism)\n")
  	cat("\t Procedure:", x$g.Ztype, "\n")
  	if(!(is.na(x$g.Zcoef[1]))){		
  		cat("\t Model:\n\t\t",x$g.Zmodel, "\n")
  		cat("\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$g.Zterms)
  		extra <- ifelse(x$g.Zcoef >= 0, "  ", " ")
  		for(i in 1:length(x$g.Zcoef)){
  			cat("\t", terms[i], extra[i], x$g.Zcoef[i], "\n")
  	}}
  	cat("\n Estimation of g.Delta (missingness mechanism)\n")
  	cat("\t Procedure:", x$g.Deltatype, "\n")
  	if(!(is.na(x$g.Deltacoef[1]))){		
  		cat("\t Model:\n\t\t",x$g.Deltamodel, "\n")
  		cat("\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$g.Deltaterms)
  		extra <- ifelse(x$g.Deltacoef >= 0, "  ", " ")
  		for(i in 1:length(x$g.Deltacoef)){
  			cat("\t", terms[i], extra[i], x$g.Deltacoef[i], "\n")
  	}}
  	cat("\n Bounds on g: (", x$gbd,")\n")
	if(!is.null(x$estimates$EY1)){
			cat("\n Population Mean")
			cat("\n   Parameter Estimate: ", signif(x$estimates$EY1$psi,5))
			cat("\n   Estimated Variance: ", signif(x$estimates$EY1$var.psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$EY1$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$EY1$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$EY1$CI[1],5), ", ", signif(x$estimates$EY1$CI[2],5), ")", sep=""),"\n") 
		}

	if(!is.null(x$estimates$ATE)){
			cat("\n Additive Effect")
			cat("\n   Parameter Estimate: ", signif(x$estimates$ATE$psi,5))
			cat("\n   Estimated Variance: ", signif(x$estimates$ATE$var.psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$ATE$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$ATE$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$ATE$CI[1],5), ", ", signif(x$estimates$ATE$CI[2],5), ")", sep=""),"\n") 
		}
		if(!is.null(x$estimates$RR)){
			cat("\n Relative Risk")
			cat("\n   Parameter Estimate: ", signif(x$estimates$RR$psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$RR$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$RR$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$RR$CI[1],5), ", ", signif(x$estimates$RR$CI[2],5), ")", sep=""),"\n") 
			cat("\n              log(RR): ", signif(x$estimates$RR$log.psi,5))
			cat("\n    variance(log(RR)): ", signif(x$estimates$RR$var.log.psi,5), "\n")

		}
		if(!is.null(x$estimates$OR)){
			cat("\n Odds Ratio")
			cat("\n   Parameter Estimate: ", signif(x$estimates$OR$psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$OR$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$OR$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$OR$CI[1],5), ", ", signif(x$estimates$OR$CI[2],5), ")", sep=""),"\n") 
			cat("\n              log(OR): ", signif(x$estimates$OR$log.psi,5))
			cat("\n    variance(log(OR)): ", signif(x$estimates$OR$var.log.psi,5),"\n")

		}
 
  } else {
  	 stop("Error calling print.summary.tmle. 'x' needs to have class 'summary.tmle'\n")
  }
}

#-------------print.tmle------------------
# print object returned by tmle function
#-----------------------------------------
print.tmle <- function(x,...) {
	if(identical(class(x), "tmle")){
		if(!is.null(x$estimates$EY1)){
			cat(" Population Mean")
			cat("\n   Parameter Estimate: ", signif(x$estimates$EY1$psi,5))
			cat("\n   Estimated Variance: ", signif(x$estimates$EY1$var.psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$EY1$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$EY1$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$EY1$CI[1],5), ", ", signif(x$estimates$EY1$CI[2],5), ")", sep=""),"\n") 
		}

		if(!is.null(x$estimates$ATE)){
			cat(" Additive Effect")
			cat("\n   Parameter Estimate: ", signif(x$estimates$ATE$psi,5))
			cat("\n   Estimated Variance: ", signif(x$estimates$ATE$var.psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$ATE$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$ATE$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$ATE$CI[1],5), ", ", signif(x$estimates$ATE$CI[2],5), ")", sep=""),"\n") 
		}
		if(!is.null(x$estimates$RR)){
			cat("\n Relative Risk")
			cat("\n   Parameter Estimate: ", signif(x$estimates$RR$psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$RR$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$RR$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$RR$CI[1],5), ", ", signif(x$estimates$RR$CI[2],5), ")", sep=""),"\n") 
			cat("\n              log(RR): ", signif(x$estimates$RR$log.psi,5))
			cat("\n    variance(log(RR)): ", signif(x$estimates$RR$var.log.psi,5),"\n")

		}
		if(!is.null(x$estimates$OR)){
			cat("\n Odds Ratio")
			cat("\n   Parameter Estimate: ", signif(x$estimates$OR$psi,5))
			cat("\n              p-value: ", ifelse(x$estimates$OR$pvalue <= 2*10^-16, "<2e-16",signif(x$estimates$OR$pvalue,5)))
			cat("\n    95% Conf Interval:",paste("(", signif(x$estimates$OR$CI[1],5), ", ", signif(x$estimates$OR$CI[2],5), ")", sep=""),"\n") 
			cat("\n              log(OR): ", signif(x$estimates$OR$log.psi,5))
			cat("\n    variance(log(OR)): ", signif(x$estimates$OR$var.log.psi,5),"\n")

		}
	} else {
 		stop("Error calling print.tmle. 'x' needs to have class 'tmle'\n")
 }}

#-------------print.tmle.list------------------
# print object returned by tmle 
# when there is a controlled direct effect
#-----------------------------------------
print.tmle.list <- function(x,...) {
	cat("Controlled Direct Effect\n")
	cat("           ----- Z = 0 -----\n")
	print(x[[1]])
	cat("\n           ----- Z = 1 -----\n")
	print(x[[2]])

}

#-------------summary.tmle.list------------------
# create tmle summary object for controlled direct 
# effect
#---------------------------------------	 
summary.tmle.list <- function(object,...) {
	summary.tmle.list <- list(Z0=summary(object[[1]]), Z1=summary(object[[2]]))
	class(summary.tmle.list) <- "summary.tmle.list"
	return(summary.tmle.list)
}

#-------------print.summary.tmle.list------------------
# create tmle summary object for controlled direct 
# effect
#---------------------------------------	 
print.summary.tmle.list <- function(x,...) {
	cat("Controlled Direct Effect\n")
	cat("           ----- Z = 0 -----\n")
	print(x[[1]])
	cat("\n           ----- Z = 1 -----\n")
	print(x[[2]])
}

# TMLE for Marginal Structural Models
# Reference: Rosenblum&vanderLaan2010, Targeted Maxium Likelihood Estimation of the Parameter of a Marginal Structural Model, The International Journal of Biostatistics, volume 6, Issue 2, Article 19.
# Note in Documentation: hAV values passed in has to be an nx2 matrix, (A=0, A=1) so that 
# we can have weights 1 for both. , same with pDelta1, at least nx2, with (A=0, A=1)
#--------
# source("/Users/Susan/cTMLE/tmle/MSM/tmleFunctionsinCommon.R")

#-------------summary.tmleMSM------------------
# create tmleMSM summary object
# return values of psi, se, pvalue, 95% conf interval,
# sigma, epsilon, psi.Qinit
# and models for Q, g, g.Delta, if available
#---------------------------------------	 
summary.tmleMSM <- function(object,...) {
  if(identical(class(object), "tmleMSM")){
		Qmodel <- Qcoef <- gmodel <- gcoef <- g.Deltamodel <- g.Deltacoef <- NULL
		g.AVmodel <- g.AVcoef <- NULL
		Qterms <- gterms <- g.Deltaterms <- g.AVterms <- ""
		if(!is.null(object$Qinit)){
		  if (!is.null(object$Qinit$coef)) {
			Qcoef <- object$Qinit$coef	
			if(class(Qcoef) == "matrix"){
				Qterms <- colnames(Qcoef)
			} else {
				Qterms <- names(Qcoef)
			}
			Qmodel <- paste("Y ~ 1")
			if(length(Qterms) > 1) {
				Qmodel <- paste("Y ~ ", paste(Qterms, collapse =" + "))
			} 
		  } else {
		  	Qmodel <- object$Qinit$type
		  }
		}
		if(!is.null(object$g)){
		if (!is.null(object$g$coef)) {
			gbd <- object$g$bound
			gcoef <- object$g$coef
			if(class(gcoef) == "matrix"){
				gterms <- colnames(gcoef)
			} else {
				gterms <- names(gcoef)
			}			
			gmodel <- paste("A ~ 1")
			if(length(gterms) > 1) {
				gmodel <- paste("A ~ ", paste(gterms, collapse =" + "))
			} 
		}}
		if(!is.null(object$g.AV)){
		if (!is.null(object$g.AV$coef)) {
			g.AVcoef <- object$g.AV$coef
			if(class(g.AVcoef) == "matrix"){
				g.AVterms <- colnames(g.AVcoef)
			} else {
				g.AVterms <- names(g.AVcoef)
			}			
			g.AVmodel <- paste("A ~ 1")
			if(length(g.AVterms) > 1) {
				g.AVmodel <- paste("A ~", paste(g.AVterms, collapse =" + "))
			} 
		}}
		if(!is.null(object$g.Delta)){
		if (!is.null(object$g.Delta$coef)) {
			g.Deltacoef <- object$g.Delta$coef
			if(class(g.Deltacoef) == "matrix"){
				g.Deltaterms <- colnames(g.Deltacoef)
			} else {
				g.Deltaterms <- names(g.Deltacoef)
			}			
			g.Deltamodel <- paste("Delta ~ 1")
			if(length(g.Deltaterms) > 1) {
				g.Deltamodel <- paste("Delta ~", paste(g.Deltaterms, collapse =" + "))
			} 
		}}
		estimates <- round(cbind(object$psi, object$se, object$pvalue, object$lb, object$ub),3)
		colnames(estimates) <- c("psi", "SE", "p-value", "lb", "ub")
		estimates[,"p-value"] <- signif(object$pvalue, 3)
		estimates[estimates[,"p-value"]< 2*10^-16, "p-value"] <- NA
		summary.tmleMSM <- list(estimates=estimates, sigma=object$sigma, 
						Qmodel=Qmodel, Qterms=Qterms, Qcoef=Qcoef, Qtype=object$Qinit$type,
						gbd=gbd, gmodel=gmodel, gterms=gterms, gcoef=gcoef,gtype=object$g$type, 
						g.AVmodel=g.AVmodel, g.AVterms=g.AVterms, g.AVcoef=g.AVcoef, g.AVtype=object$g.AV$type,
						g.Deltamodel=g.Deltamodel, g.Deltaterms=g.Deltaterms, g.Deltacoef=g.Deltacoef, g.Deltatype=object$g.Delta$type, psi.Qinit=object$psi.Qinit)
		class(summary.tmleMSM) <- "summary.tmleMSM"
	} else {
		stop("object must have class 'tmleMSM'")
		summary.tmle <- NULL
	}
	return(summary.tmleMSM)
}

#-------------print.summary.tmleMSM------------------
# print tmleMSM summary object
#-------------------------------------------------
print.summary.tmleMSM <- function(x,...) {
  if(identical(class(x), "summary.tmleMSM")){
  	cat(" Initial estimation of Q\n")
  	cat("\t Procedure:", x$Qtype)
  	if(!(is.na(x$Qcoef[1]))){
  		cat("\n\t Model:\n\t\t", x$Qmodel)
  		cat("\n\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$Qterms)
  		extra <- ifelse(x$Qcoef >= 0, "  ", " ")
  		for(i in 1:length(x$Qcoef)){
  			cat("\t", terms[i], extra[i], x$Qcoef[i], "\n")
  		}
  	}
  	cat("\n Estimation of g (treatment mechanism)\n")
  	cat("\t Procedure:", x$gtype,"\n")
  	if(!(is.na(x$gcoef[1]))){	
 		cat("\t Model:\n\t\t", x$gmodel, "\n")
  		cat("\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$gterms)
  		extra <- ifelse(x$gcoef >= 0, "  ", " ")
  		for(i in 1:length(x$gcoef)){
  			cat("\t", terms[i], extra[i], x$gcoef[i], "\n")
  	}}
  	cat("\n Estimation of h(A,V)\n")
  	cat("\t Procedure:", x$g.AVtype, "\n")
  	if(!(is.na(x$g.AVcoef[1]))){		
  		cat("\t Model:\n\t\t",x$g.AVmodel, "\n")
  		cat("\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$g.AVterms)
  		extra <- ifelse(x$g.AVcoef >= 0, "  ", " ")
  		for(i in 1:length(x$g.AVcoef)){
  			cat("\t", terms[i], extra[i], x$g.AVcoef[i], "\n")
  	}}
  	cat("\n Estimation of g.Delta (missingness mechanism)\n")
  	cat("\t Procedure:", x$g.Deltatype, "\n")
  	if(!(is.na(x$g.Deltacoef[1]))){		
  		cat("\t Model:\n\t\t",x$g.Deltamodel, "\n")
  		cat("\n\t Coefficients: \n")
  		terms <- sprintf("%15s", x$g.Deltaterms)
  		extra <- ifelse(x$g.Deltacoef >= 0, "  ", " ")
  		for(i in 1:length(x$g.Deltacoef)){
  			cat("\t", terms[i], extra[i], x$g.Deltacoef[i], "\n")
  	}}
  	cat("\n Bounds on g: (", x$gbd,")\n")
	if(!is.null(x$estimates)){
		cat("\n MSM parameter estimates and 95% confidence intervals\n")
		print(x$estimates,quote=FALSE, na.print="<2e-16")
	}
	if(!is.null(x$sigma)){
		cat("\n Variance-covariance matrix\n")
		print(x$sigma,3)
	} 
  } else {
  	 stop("Error calling print.summary.tmleMSM. 'x' needs to have class 'summary.tmleMSM'\n")
  }
}

#-------------print.tmleMSM------------------
# print object returned by tmleMSM function
#-----------------------------------------
print.tmleMSM <- function(x,...) {
	if(identical(class(x), "tmleMSM")){
		if(!is.null(x$psi) & !is.null(x$sigma) & !is.null(x$lb)){
			cat(" MSM parameter estimates and 95% confidence intervals\n")
			estimates <- round(cbind(x$psi, x$se, x$pvalue, x$lb, x$ub),3)
			colnames(estimates) <- c("psi", "SE", "p-value", "lb", "ub")
			estimates[,"p-value"] <- signif(x$pvalue, 3)
			estimates[estimates[,"p-value"]< 2*10^-16, "p-value"] <- NA
			print(estimates,3,quote=FALSE, na.print="<2e-16")
		} else {
			if(!is.null(x$psi)){
		    	cat(" MSM parameters\n  ")
		    	print(x$psi,3)
			}
			if(!is.null(x$sigma)){
		    	cat("\n Standard Errors\n  ")
		    	names(x$se) <- colnames(x$sigma)
		    	print(x$se,3)
			} 
			if(!is.null(x$lb)){
		    	cat("\n 95% Confidence Intervals\n  ")
		    	print(rbind(lb=x$lb, ub=x$ub),3)
			} 
		}	
	} else {
 		stop("Error calling print.tmle. 'x' needs to have class 'tmle'\n")
 }}

# This is defined differently in SL v 2.0+, undefined in SL <2.0
SL.glm.interaction <- function (Y.temp, X.temp, newX.temp, family, obsWeights, ...) 
{
    fit.glm <- glm(Y.temp ~ .^2, data = X.temp, family = family, weights = obsWeights)
    out <- predict(fit.glm, newdata = newX.temp, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    foo <- list(out=out, fit = fit)
    return(foo)
}


# Computes the variance/covariance matrix for all params in the IC.
calcSigma <- function(hAV, gAVW, Y, Q, mAV, covar.MSM, covar.MSMA0, covar.MSMA1, I.V, Delta, ub, id, family){
	D <- I.V * Delta* ( .bound(hAV[,"hAV"]/gAVW, c(0, ub)) * (Y-Q[,"QAW"]) * covar.MSM
		+ hAV[,"h0V"]* (Q[,"Q0W"] - mAV[,"m0V"]) * covar.MSMA0
		+ hAV[,"h1V"]* (Q[,"Q1W"] - mAV[,"m1V"]) * covar.MSMA1)
	if(!any(is.na(D))){
   		# Construct normalizing matrix M
   		nterms <- ncol(covar.MSMA0)
   		f <- function(x){x[1:nterms] %*% t(x[(nterms+1):(2*nterms)])}
   		if(family == "binomial"){
   			derivFactor <- cbind(mAV[,"m0V"] * (1-mAV[,"m0V"]), mAV[,"m1V"] * (1-mAV[,"m1V"]))
   		} else {
   			derivFactor <- matrix(1, nrow=nrow(mAV), ncol = 2)
   		}
   		deriv.term2 <-  apply(cbind(-hAV[,"h0V"]* I.V * derivFactor[,1] * covar.MSMA0, covar.MSMA0), 1, f)
   		deriv.term3 <-  apply(cbind(-hAV[,"h1V"]* I.V * derivFactor[,2] * covar.MSMA1, covar.MSMA1), 1, f) 
   		ddpsi.D <- as.matrix(deriv.term2 + deriv.term3)
   		M <- -matrix(rowMeans(ddpsi.D), nrow=nterms)

   		Minv <- try(solve(M))
   		if(identical(class(Minv), "try-error")){
   			warning("Inference unavailable: normalizing matrix not invertible\n")
   			sigma <- NA
   		} else {
   			Dstar <- t(Minv %*% t(D))
   			if(length(unique(id)) < length(Y)){
   				Dstar <- matrix(unlist(by(Dstar, id, colMeans, simplify=TRUE)), byrow=TRUE, nrow=length(unique(id)))
			}
   			sigma <- var(Dstar)
   		}
   		rownames(sigma) <- colnames(sigma) <- colnames(covar.MSMA0)
   	} else {
   		D <- .bound(hAV[,"hAV"]/gAVW, c(0, ub)) * (Y-Q[,"QAW"]) * covar.MSM
		term2 <-  hAV[,"h0V"]* (Q[,"Q0W"] - mAV[,"m0V"]) * covar.MSMA0
		term3 <- hAV[,"h1V"]* (Q[,"Q1W"] - mAV[,"m1V"]) * covar.MSMA1
		if (!any(is.na(term2))){D <- D + term2}
		if (!any(is.na(term3))){D <- D + term3}
		sigma <- var(I.V * Delta * D)
	}
   return(sigma)
}

tmleMSM <- function(Y,A,W,V,T=rep(1,length(Y)), Delta=rep(1, length(Y)), MSM, v=NULL, 
				Q=NULL, Qform=NULL, 
				Qbounds=c(-Inf, Inf), Q.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"), cvQinit = FALSE,
				hAV=NULL, hAVform=NULL,  
				g1W = NULL, gform=NULL, 
				pDelta1=NULL, g.Deltaform=NULL, 
				g.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"), ub = 1/0.025, 
				family="gaussian", fluctuation="logistic", alpha  = 0.995,
				id=1:length(Y), inference=TRUE, verbose=FALSE) {
	Y[is.na(Y)] <- 0
	n <- length(Y)
	n.id <- length(unique(id))
	if(is.null(v)){
		I.V <- rep(1, length(Y))
	} else {
		I.V <- as.numeric(V==v)
	}
	V <- as.matrix(V)
	colnames(V) <- .setColnames(colnames(V), NCOL(V), "V")

	if (sum(!sapply(W, is.numeric)) > 0) {
        stop("Currently, only numeric variables are allowed.  Please convert any character or factor variables to numeric.")
    }

	W <- as.matrix(W)
	colnames(W) <- .setColnames(colnames(W), NCOL(W), "W")
	
	if(identical (family, binomial)){
		family <- "binomial"
	} else if (identical(family, gaussian)){
		family <- "gaussian"
	}
	if(!.verifyArgs(Y,Z=NULL,A,cbind(V,W,T),Delta, Qform, gform, hAVform, g.Deltaform)){
		stop()
	}
	maptoYstar <- fluctuation=="logistic" | family=="binomial"

#---- Stage 1 -----
	stage1 <- .initStage1(Y=Y,A=A, Q=Q, Delta=Delta, Qbounds=Qbounds, alpha=alpha, maptoYstar=maptoYstar, family=family)
	Qinit <- suppressWarnings(estimateQ(Y=stage1$Ystar,Z=rep(1, length(Y)), A=A, 
			W=cbind(W,V,T), Delta=(I.V==1 & Delta==1),
			Q=stage1$Q, Qbounds=stage1$Qbounds, Qform=Qform, maptoYstar = maptoYstar,
			SL.library=Q.SL.library, cvQinit=cvQinit, family=family, id=id, verbose=verbose))

#---- Stage 2 -----
    if(is.null(hAV)){
    	gAV <- suppressWarnings(estimateG(d=data.frame(A,V,T), hAV, hAVform,g.SL.library, id, verbose, 
    			message="h(A,V)", outcome="A"))
		hAV <- cbind((1-A)*(1-gAV$g1W) + A*gAV$g1W, 1-gAV$g1W, gAV$g1W)
	} else {
		hAV <- cbind((1-A)*(hAV[,1]) + A*hAV[,2], hAV)
		gAV <- NULL
		gAV$g1W <- hAV
		gAV$type <- "User-supplied values"
		gAV$coef <- NA
	}
	colnames(hAV) <- c("hAV", "h0V", "h1V")
	if (is.null(v)){
		g <- suppressWarnings(estimateG(d=data.frame(A,V,W,T), g1W, gform,g.SL.library, id, verbose, 
    			message="treatment mechanism", outcome="A"))
	} else {	
		g <- suppressWarnings(estimateG(d=data.frame(A,V,W,T), g1W, gform,g.SL.library, id, verbose, 
    			message="treatment mechanism", outcome="A", newdata=data.frame(A,V=v, W,T)))
   }
   g$bound <- c(0,ub)
   if(g$type=="try-error"){
 		stop("Error estimating treatment mechanism (hint: only numeric variables are allowed)") 
 	}
 	g.Delta <- estimateG(d=data.frame(Delta, Z=1, A, W,V,T), pDelta1, g.Deltaform, 
 		g.SL.library,id=id, verbose, "missingness mechanism", outcome="D") 
	g1VW <- g$g1W * g.Delta$g1W[,"Z0A1"]
	g0VW <- (1-g$g1W) * g.Delta$g1W[,"Z0A0"]
	gAVW <- A*g1VW + (1-A)*g0VW

	MSMformula <- formula(paste("Y~",MSM))
	
	mfA <- model.frame(MSMformula, data=data.frame(Y,A,V,W,T))  # Added Y Oct 24 ,2011
	mf0 <- model.frame(MSMformula, data=data.frame(Y,A=rep(0,n),V,W,T))
	mf1 <- model.frame(MSMformula, data=data.frame(Y,A=rep(1,n),V,W,T))
	
	covar.MSM   <- model.matrix(MSMformula, mfA)
	covar.MSMA0 <- model.matrix(MSMformula, mf0)
	covar.MSMA1 <- model.matrix(MSMformula, mf1)
		
	if(verbose){cat("\tTargeting Q\n")}
	C1 <- I.V * .bound(hAV[,"hAV"]/gAVW, c(0,ub)) * covar.MSM  
	suppressWarnings(
	  epsilon <- coef(glm(stage1$Ystar ~ -1 + offset(Qinit$Q[,"QAW"]) + C1, subset=(I.V==1 & Delta==1), family=Qinit$family))
	)
	Qstar <- cbind(Qinit$Q[,"QAW"] + C1 %*% epsilon,
					 Qinit$Q[,"Q0W"] + .bound(hAV[,"h0V"]/g0VW, c(0,ub)) * covar.MSMA0 %*% epsilon,
					 Qinit$Q[,"Q1W"] + .bound(hAV[,"h1V"]/g1VW, c(0,ub)) * covar.MSMA1 %*% epsilon)
	if(identical(Qinit$family, "binomial")){
		Qstar <- plogis(Qstar)*diff(stage1$ab)+stage1$ab[1] 
		Qinit$Q <- plogis(Qinit$Q)*diff(stage1$ab)+stage1$ab[1]
		Ystar <- stage1$Ystar*diff(stage1$ab)+stage1$ab[1]
	}
   colnames(Qstar) <- c("QAW", "Q0W", "Q1W")
  
   if(verbose){cat("\tEvaluating MSM parameters\n")}
	d.Qstar <- data.frame(Y=c(Qstar[,"Q0W"], Qstar[,"Q1W"]), 
							rbind(mf0, mf1),
							wts=c(hAV[,"h0V"], hAV[,"h1V"]))
	suppressWarnings(
		psi.Qstar <- coef(glm(MSMformula, data=d.Qstar, 
							weights=d.Qstar$wts, family=family)) 
	)
	d.Qinit <- replace(d.Qstar,1, c(Qinit$Q[,"Q0W"], Qinit$Q[,"Q1W"]))
	suppressWarnings(
		psi.Qinit <- coef(glm(MSMformula, data=d.Qinit, weights=d.Qinit$wts, family=family))
	)

    if(inference){
    	if(verbose){cat("\tCalculating variance-covariance matrix\n")}
    	if(family=="binomial"){
    		mAV <- plogis(cbind(covar.MSMA0 %*% psi.Qstar, covar.MSMA1 %*% psi.Qstar))
    	} else {
    		mAV <- matrix(1, nrow=nrow(covar.MSMA0), ncol=2)
    	}
    	colnames(mAV) <- c("m0V", "m1V")
    	sigma <- calcSigma(hAV, gAVW, Y, Qstar, mAV, covar.MSM, covar.MSMA0, covar.MSMA1, I.V, Delta, ub, id, family)/n.id
    	se <- sqrt(diag(sigma))
	pvalue <- 2*pnorm(-abs(psi.Qstar/se))
    	lb <- psi.Qstar -1.96*se
    	ub <- psi.Qstar +1.96*se
    } else {
    	sigma <- se <- lb <- ub <- NULL
    }
    Qinit$Q <- Qinit$Q[,-1]
    returnVal <- list(psi=psi.Qstar, sigma=sigma,se=se, pvalue=pvalue, lb=lb, ub=ub, epsilon=epsilon,  psi.Qinit=psi.Qinit,  Qstar=Qstar[,-1], Qinit=Qinit, g=g, g.AV=gAV, g.Delta=g.Delta)
    class(returnVal) <- "tmleMSM" 
	return(returnVal)
}

#-------------.verifyArgs------------------
# initial checks on data passed in
#-------------------------------------------
.verifyArgs <- function(Y,Z,A,W,Delta, Qform, gform, g.Zform,g.Deltaform){
	formulas <- list(Qform, gform, g.Zform, g.Deltaform)
	validFormula <- sapply(formulas, function(x){identical(class(try(as.formula(x))), "formula")})
	validNames <- c("Y", "Z", "A", ".", "Delta", colnames(W))
	validTerms <- rep(TRUE, length(formulas))
	validTerms[validFormula] <- sapply(formulas[which(validFormula)], 
			function(x){is.null(x) || all(all.names(as.formula(x), functions=FALSE) %in% validNames)})
	ok <- c(length(Y) == length(A) & length(A) == NROW(W) & length(A)==NROW(Delta),
			sum(is.na(A), is.na(W)) == 0,
	 	 	all(A[!is.na(A)] %in% 0:1),
	 	 	is.null(Z) || all(Z[!is.na(Z)] %in% 0:1),
	 	 	is.null(Z) || length(unique(Z)) == 1 | (length(unique(Z)) > 1 & length(unique(A))>1),
			validFormula,
			validTerms
			)
	warning_messages <- c("\t'Y', 'A', 'W', 'Delta', must contain the same number of observations\n",
				"\tNo missing values allowed in 'A' or 'W'\n",
				"\t'A' must be binary (0,1)\n",
				"\t'Z' must be binary (0,1)\n",
				"\tIntermediate variable (Z) not allowed when there is no experimentation in A",
				"\tInvalid regression formula for 'Qform'",
				"\tInvalid regression formula for 'gform'",
				"\tInvalid regression formula for 'g.Zform'",
				"\tInvalid regression formula for 'g.Deltaform'",
				"\tInvalid term name in regression formula for 'Qform'",
				"\tInvalid term name in regression formula for 'gform'",
				"\tInvalid term name in regression formula for 'g.Zform'",
				"\tInvalid term name in regression formula for 'g.Deltaform'"
				)
	if(!all(ok)){
		warning("\n", warning_messages[!ok], immediate. = TRUE)
	}
	return(all(ok))
}

#---------- function .setColnames ---------------
# assign names to every unnamed column of x
# arguments
# 	x.colnames - current column names
#	x.ncols - current number of columns
# 	firstChar - prefix for internally assigned name
# return the names
#-----------------------------------------
.setColnames <- function(x.colnames, x.ncols, firstChar){
	if(is.null(x.colnames)) {
		if(x.ncols > 1){
			x.colnames <- paste(firstChar,1:x.ncols, sep="")
		} else {
			x.colnames <- firstChar
		}
	} else {
		invalid.name <- nchar(x.colnames) == 0
		if(any(invalid.name)){
			x.colnames[invalid.name] <- paste(".internal",firstChar, which(invalid.name), sep="")
		}
	}
	return(x.colnames)
}

#---------- function .bound ---------------
# set outliers to min/max allowable values
# assumes x contains only numerical data
#-----------------------------------------
.bound <- function(x, bounds){
	x[x>max(bounds)] <- max(bounds)
	x[x<min(bounds)] <- min(bounds)
	return(x)
}

#---------- function .initStage1 ---------------
# Bound Y, map to Ystar if applicable, and
# set boundson on Q and enforce on user-specified values
# returns
#   Ystar - outcome values (between [0,1] if maptoYstar=TRUE)
#   Q - matrix of user-specified values
#   Qbounds - bounds on predicted values for Q (10% wider at each end then
# 			observed range of Y
#			(-Inf,+Inf) is default for linear regression
#   ab - bounding levels used to transform Y to Ystar
#-----------------------------------------------
.initStage1 <- function(Y,A, Q, Q.Z1=NULL, Delta, Qbounds, alpha, maptoYstar, family){
	if(family=="binomial") {Qbounds <- c(0,1)}
 	if(is.null(Qbounds)) {
 		if(maptoYstar){ 
 			Qbounds <- range(Y[Delta==1])
 			Qbounds <- Qbounds + .1*c(-abs(Qbounds[1]),abs(Qbounds[2]))
 		} else {
 			Qbounds <- c(-Inf, Inf)
 		}
   	}
	if(!is.null(Q)){
   		QAW <- (1-A)*Q[,1] + A*Q[,2]
   		Q <- cbind(QAW, Q0W=Q[,1], Q1W=Q[,2])
   	} 
   	if(!is.null(Q.Z1)){
   		Q <- cbind(Q, Q0W.Z1=Q.Z1[,1], Q1W.Z1=Q.Z1[,2])
   	}
  	ab <- c(0,1)
   	Ystar <- Y
   	if(maptoYstar){ 
   		Ystar <- .bound(Y, Qbounds)
   		if(!is.null(Q)){
   			Q <- .bound(Q, Qbounds)
   		}
   		if(0 >= alpha | 1 <= alpha){
			alpha <- .995
			warning(paste("\n\talpha must be between 0 and 1, alpha reset to",alpha,"\n"),
							immediate. = TRUE)
		}
		ab <- range(Ystar, na.rm=TRUE)
		Ystar[is.na(Ystar)] <- 0  
		Ystar <- (Ystar-ab[1])/diff(ab)	
		if(!is.null(Q)){Q <- (Q-ab[1])/diff(ab)}
		Qbounds <- c(alpha, 1-alpha)
	}	
	return(list(Ystar=Ystar, Q=Q, Qbounds=Qbounds, ab=ab))
} 


#----- function .estQcvSL ----
# purpose: Obtain cross-validated estimates for initial Q using discrete SL.  
# 	This function will return cross-validated predicted initial values 
#   corresponding to the best algorithm in the library, or the convex combination
#   calculated by SL itself, as chosen by cvRSS.
#   The fitted value for each observation i, is predicted from a fit based on a training set
#	excluding the fold containing observation i. Observations with same id are grouped in the
#   same fold.
# arguments: 
# 	Y - outcome
# 	X - design matrix with (Z,A,W)
# 	SL.library - prediction algorithms
# 	V - number of outer cv-folds for disscrete SL
# 	V_SL - number of folds for internal SL cross-validation
# 	family - binomial or gaussian
#   Delta - missingness indicator
#   Qbounds - bounds on predicted values for Q
#   id - subject identifier
# returns:
#  Q - nx5 matrix of predicted values on linear scale (e.g. logit if Qfamily=binomial)
#-----------------------------------------------
.estQcvSL <- function(Y,X,SL.library=NULL, V=5, V_SL=5, family="gaussian", Delta, Qbounds, id, verbose){
	SL.version <- 2
	Q <- cvRSS <- best_alg <- NULL
	n <- length(Y)
	u.id <- unique(id)
	n.id <- length(u.id)
	fold <- by(sample(1:n.id),rep(1:V, length.out=n.id),function(x){which(id %in% u.id[x])})
	n_predictors <-length(SL.library)
	CDE <- length(unique(X[,1])) > 1

	if(NCOL(X) > 1){
		newX <- rbind(X,X,X,X,X)
		newX[(n+1):(3*n),1]   <- 0 
		newX[(3*n+1):(5*n),1] <- 1
		newX[(n+1):(2*n),2]   <- 0
		newX[(2*n+1):(3*n),2] <- 1
		newX[(3*n+1):(4*n),2] <- 0
		newX[(4*n+1):(5*n),2] <- 1
		# We'll create a matrix of predictions - one column for each predictor in the library
		# plus one more for SL itself, with 5*n predicted values per column, corresponding to newX.
   		predictions <- matrix(nrow=5*n, ncol=length(SL.library)+1)
   		m_SL <- NULL
    	for (v in 1:V) {
    		fold_rows <- c(fold[[v]], fold[[v]]+n, fold[[v]]+2*n, fold[[v]]+3*n, fold[[v]]+4*n)
    		if(class(m_SL) != "try-error"){
    			train.observed <- (1:n)[-fold[[v]]][Delta[-fold[[v]]]==1]
    			if(packageDescription("SuperLearner")$Version < SL.version){
    				arglist <- list(Y=Y[train.observed], X=X[train.observed,], newX=newX[fold_rows,],
    			 	V=V_SL, save.fit.library=FALSE, family=family,SL.library=SL.library,id=id[train.observed])
    			} else {
    				arglist <- list(Y=Y[train.observed], X=X[train.observed,], newX=newX[fold_rows,],
    			 	cvControl=list(V=V_SL), control = list(saveFitLibrary=FALSE), family=family,SL.library=SL.library,id=id[train.observed])
    			}

    			suppressWarnings(
    				m_SL <- try(do.call(SuperLearner, arglist))
    			 )
    		}
    		if(class(m_SL) != "try-error"){
    			predictions[fold_rows,1] <- m_SL$SL.predict
    			for (s in 1:n_predictors){
    				predictions[fold_rows,s+1] <- m_SL$library.predict[,s]
    			}
    			predictions <- .bound(predictions, Qbounds)
    		} else {
    			stop("Super Learner failed when estimating Q. Exiting program\n")
    		}
    	} 
    	cvRSS <- colSums(Delta*(Y-predictions[1:n,])^2)
    	names(cvRSS) <- c("SL", SL.library)
    	best <- which.min(cvRSS)
    	best_alg <- c("SL", SL.library)[best]
    	Q <- matrix(data=predictions[,best], nrow=n, ncol=5, byrow=FALSE)
    	colnames(Q) <- c("QAW", "Q0W", "Q1W", "Q0W.Z1", "Q1W.Z1")
    } 	
	if(verbose){cat("\tDiscrete SL: best algorithm = ", best_alg,"\n")}
	if (is.null(Q) | class(m_SL) == "try-error"){
		Q <- 0
		class(Q) <- "try-error"
	}
	Qinit <- list(Q=Q, family=family, SL.library=SL.library, cvRSS=cvRSS, best_alg=best_alg)
    return(Qinit)
 }

#-----------estimateQ----------------
# purpose: estimate Q=E(Y |Z, A,W) data-adaptively,
# unless super learner not available, or user specifies 
# initial values or a regression formula
# arguments: 
# 	Y - outcome 
# 	Z - intermediate variable between A and Y (default= 0 when no int. var.) 
#	A - treatment indicator (1=treatment, 0=control)
# 	W - baseline covariates
#	Delta - missingness indicator
#	Q - optional externally estimated values for Q
#	Qbounds - bounds for predicted values 
#  	Qform - optional regression formula to use for glm if 
#	        non-data adaptive estimation specified
# 	maptoYstar - if TRUE, using logistic fluctuation for bounded, continuous outcomes
# 		estimation inital Q on linear scale, bounded by (0,1),and return on logit scale
#		(will work if family=poisson)
#	SL.library - library of prediction algorithms for Super Learner
#   cvQinit - flag, if TRUE, cross-validate SL.
# 	family - regression family
#	id - subject identifier
# returns matrix of linear predictors for Q(A,W), Q(0,W), Q(1,W),
#   (for controlled direct effect, 2 additional columns: Q(Z=1,A=0,W), Q(Z=1,A=1,W)) 
#		family for stage 2 targeting
#		coef, NA, unless Q is estimated using a parametric model
# 		type, estimation method for Q
#----------------------------------------
estimateQ <- function (Y,Z,A,W, Delta, Q, Qbounds, Qform, maptoYstar, 
		SL.library, cvQinit, family, id, verbose) {
	SL.version <- 2
	Qfamily <- family
	m <- NULL
	coef <- NA
	CDE <- length(unique(Z)) > 1
	type <- "user-supplied values"
	if(is.null(Q)){
		if(verbose) { cat("\tEstimating initial regression of Y on A and W\n")}
		Q <- matrix(NA, nrow=length(Y), ncol = 5)
  	  	colnames(Q)<- c("QAW", "Q0W", "Q1W", "Q0W.Z1", "Q1W.Z1")
  	  	if(!(is.null(Qform))){
  	  		if(identical(as.character(as.formula(Qform)), c("~","Y", "."))){
  	  			if(CDE){
  	  				Qform <- paste("Y~Z+A+", paste(colnames(W), collapse="+"))
  	  			} else {
  	   	  			Qform <- paste("Y~A+", paste(colnames(W), collapse="+"))
				}
  	  		}
  	  		m <- suppressWarnings(glm(Qform, data=data.frame(Y,Z,A,W, Delta), family=family, subset=Delta==1))
  	  		Q[,"QAW"] <- predict(m, newdata=data.frame(Y,Z,A,W), type="response")
  	  		Q[,"Q0W"] <- predict(m, newdata=data.frame(Y,Z=0,A=0,W), type="response")
  	  		Q[,"Q1W"] <- predict(m, newdata=data.frame(Y,Z=0,A=1,W), type="response")
  	  		Q[,"Q0W.Z1"] <- predict(m, newdata=data.frame(Y,Z=1,A=0,W), type="response")
  	  		Q[,"Q1W.Z1"] <- predict(m, newdata=data.frame(Y,Z=1,A=1,W), type="response")
  	  		coef <- coef(m)
  	  		type="glm, user-supplied model"	
  	  	} else {
  	  		if(cvQinit){
  	  			m <- try(.estQcvSL(Y,X=cbind(Z,A,W),SL.library, family=family, 
  	  					Delta=Delta, Qbounds=Qbounds,id=id, verbose=verbose))
  	  			if(!(identical(class(m), "try-error"))){
  	  				type <- "cross-validated SL"
  	  				Qinit <- m
  	  				Q <- Qinit$Q
  	  			}
  	  		} else {
  				if(verbose) {cat("\t using SuperLearner\n")}
  				n <- length(Y)
  				X <- data.frame(Z,A,W)
  				X00 <- data.frame(Z=0,A=0, W)
  				X01 <- data.frame(Z=0,A=1, W)
  				newX <- rbind(X, X00, X01)
  				if(CDE) {
  					X10 <- data.frame(Z=1,A=0, W)
  					X11 <- data.frame(Z=1,A=1, W)
  					newX <- rbind(newX, X10, X11)
  				}
  				if(packageDescription("SuperLearner")$Version < SL.version){
    				arglist <- list(Y=Y[Delta==1],X=X[Delta==1,], newX=newX, SL.library=SL.library,
  							V=5, family=family, save.fit.library=FALSE, id=id[Delta==1])
  				} else {
    				arglist <- list(Y=Y[Delta==1],X=X[Delta==1,], newX=newX, SL.library=SL.library,
    			 		cvControl=list(V=5), family=family, control = list(saveFitLibrary=FALSE), id=id[Delta==1])
    				}
    				suppressWarnings(
    					m<- try(do.call(SuperLearner, arglist))
    			 	)	
  				if(identical(class(m),"SuperLearner")){
  					Q[,"QAW"] <- m$SL.predict[1:n]
  					Q[,"Q0W"] <- m$SL.predict[(n+1):(2*n)]
  					Q[,"Q1W"] <- m$SL.predict[(2*n+1):(3*n)]
  					if(CDE){
  						Q[,"Q0W.Z1"] <- m$SL.predict[(3*n+1):(4*n)]
  				  		Q[,"Q1W.Z1"] <- m$SL.predict[(4*n+1):(5*n)]
  					}
  					type <- "SuperLearner"
  	  			} else {
  	  				stop("Super Learner failed when estimating Q. Exiting program\n")
  	  			} 
  	   } }
  	} 
  	if(is.na(Q[1,1]) | identical(class(m), "try-error")){
  	  		if(verbose) {cat("\t Running main terms regression for 'Q' using glm\n")}
  	  		Qform <- paste("Y~Z+A+", paste(colnames(W), collapse="+"))
  	  		m <- glm(Qform, data=data.frame(Y,Z,A,W, Delta), family=family, subset=Delta==1)
  	  		Q[,"QAW"] <- predict(m, newdata=data.frame(Y,Z,A,W), type="response")
  	  		Q[,"Q1W"] <- predict(m, newdata=data.frame(Y,Z=0,A=1,W), type="response")
  	  		Q[,"Q0W"] <- predict(m, newdata=data.frame(Y,Z=0,A=0,W), type="response")
	  		Q[,"Q0W.Z1"] <- predict(m, newdata=data.frame(Y,Z=1,A=0,W), type="response")
  	  		Q[,"Q1W.Z1"] <- predict(m, newdata=data.frame(Y,Z=1,A=1,W), type="response")

  	  		coef <- coef(m)
  	  		type="glm, main terms model"
  	 }
	Q <- .bound(Q, Qbounds)
	if(maptoYstar | identical(Qfamily,"binomial") | identical(Qfamily, binomial)){
			Q <- qlogis(Q)
			Qfamily <- "binomial"	
	} else if (identical(Qfamily, "poisson") | identical(Qfamily, poisson)) {
			Q <- log(Q)
			Qfamily <- "poisson"
	}
	if(!CDE){
		Q <- Q[,1:3]
	}
	if(cvQinit){
		Qinit$Q <- Q
	} else {
		Qinit <- list(Q=Q, family=Qfamily, coef=coef, type=type)
		if(type=="SuperLearner"){
			 Qinit$SL.library=SL.library
			 Qinit$coef=m$coef
		}
	}
	return(Qinit)
}

#-----------estimateG----------------
# Estimate factors of g
# 		P(A=1|W), P(Z=1|A,W), P(Delta=1|Z,A,W)
# d - dataframe (A,W), (Z,A,W), or (Delta,Z,A,W)
# g1W - optional vector/matrix  of externally estimated values
# gform - optionalformula to use for glm
# SL.library - algorithms to use for super learner estimation
#   id - subject identifier
# verbose - flag, whether or not to print messages
# message - printed when verbose=TRUE 
# outcome - "A" for treatment, "Z" for intermediate variable,
#           "D" for Delta (missingness)
# newdata - optional values to predict on (needed by tmleMSM function)
# d = [A,W] for treatment
# d = [Z,A,W] for intermediate
# d = [Delta, Z,A,W for missingness]
#----------------------------------------
estimateG <- function (d,g1W, gform,SL.library, id, verbose, message, outcome="A", newdata=d)  {
  SL.version <- 2
  SL.ok <- FALSE
  m <- NULL
  coef <- NA
  type <- NULL
  if (is.null(g1W)){
  	if(verbose){cat("\tEstimating", message, "\n")}
	if (length(unique(d[,1]))==1) {
		g1W <- rep(1,nrow(d))
		type <- paste("No", strsplit(message, " ")[[1]][1])
		if(outcome=="Z"){
			g1W <- cbind(A0=g1W, A1=g1W)
		} else if (outcome=="D"){
			g1W <- cbind(Z0A0=g1W, Z0A1=g1W, Z1A0=g1W, Z1A1=g1W)
		}
	} else {
  	  if (is.null(gform)){
  	  	SL.ok <- TRUE
  		old.SL <- packageDescription("SuperLearner")$Version < SL.version
  		if(old.SL){
  			arglist <- list(Y=d[,1], X=d[,-1, drop=FALSE], newX=newdata[,-1, drop=FALSE], family="binomial", SL.library=SL.library, V=5, id=id)
  		} else {
  		 	arglist <- list(Y=d[,1], X=d[,-1, drop=FALSE], newX=newdata[,-1, drop=FALSE], family="binomial", SL.library=SL.library, cvControl=list(V=5), id=id)
  		}
  		suppressWarnings(
  			m <- try(do.call(SuperLearner,arglist))
  		)
  		if(identical(class(m),"SuperLearner")) {
  			g1W <- as.vector(m$SL.predict)
  		} else {
  			SL.ok <- FALSE
  			cat("Error estimating g using SuperLearner. Defaulting to glm\n")
  		}
  	    if (!SL.ok){
  			if(verbose){cat("\tRunning main terms regression for 'g' using glm\n")}
			form <- paste(paste(colnames(d)[1],"~1"), paste(colnames(d)[-1], collapse = "+"), sep="+")  
			m <- glm(form, data=d, family="binomial")
			g1W <- predict(m, newdata=newdata, type="response")
			coef <- coef(m)
  		}
  	  } else {
  	  	form <- try(as.formula(gform))
  	  	if(class(form)== "formula") {
  	  		m <- try(glm(form,  data=d, family="binomial"))
  	  		if (class(m)[1]=="try-error"){
  				if(verbose){cat("\tInvalid formula supplied. Running glm using main terms\n")}
				form <- paste(colnames(d)[1],"~1 + ", paste(colnames(d)[-1], collapse = "+"), sep="")  
  	  			m <- glm(form, data=d, family="binomial")
  	  		} else {
  	  			type <- "user-supplied regression formula"
  	  		}
  	  	} else {
  	  	  	if(verbose){cat("\tRunning main terms regression for 'g' using glm\n")}
			form <- paste(colnames(d)[1],"~1", paste(colnames(d)[-1], collapse = "+"), sep="+")  
  	  		m <- glm(form, data=d, family="binomial")
  	  	}
  	  	 g1W <- predict(m, newdata=newdata, type="response")
  	  	 coef <- coef(m)
  	  }
  	  # Get counterfactual predicted values
  	  if(outcome=="Z"){
  	  	if(identical(class(m),"SuperLearner")){
  	  		g1W <- cbind(predict(m, newdata=data.frame(A=0, newdata[,-(1:2), drop=FALSE]), type="response", 
  	  								X=d[,-1, drop=FALSE], Y=d[,1])[[1]], 
  	  				      predict(m, newdata=data.frame(A=1, newdata[,-(1:2), drop=FALSE]), type="response", 
  	  				   				X=d[,-1, drop=FALSE], Y=newdata[,1])[[1]])
  	  	} else {
  	  		g1W <- cbind(predict(m, newdata=data.frame(A=0, newdata[,-(1:2), drop=FALSE]), type="response"), 
  	  				   predict(m, newdata=data.frame(A=1, newdata[,-(1:2), drop=FALSE]), type="response"))
  	  	}	
  	  	colnames(g1W) <- c("A0", "A1")
		
  	  } else if (outcome=="D"){
  	  	if(identical(class(m),"SuperLearner")){
  	  		g1W <- cbind(predict(m, newdata=data.frame(Z=0, A=0, newdata[,-(1:3), drop=FALSE]), type="response", 
  	  							X=d[,-1,drop=FALSE], Y=d[,1])[[1]], 
  	  			 		  predict(m, newdata=data.frame(Z=0, A=1, newdata[,-(1:3), drop=FALSE]), type="response", 
  	  			 				X=d[,-1, drop=FALSE], Y=d[,1])[[1]],
  	  		  	 		  predict(m, newdata=data.frame(Z=1, A=0, newdata[,-(1:3), drop=FALSE]), type="response", 
  	  		  	 				X=d[,-1, drop=FALSE], Y=d[,1])[[1]],	     	     		  
  	  		  	 			predict(m, newdata=data.frame(Z=1, A=1, newdata[,-(1:3), drop=FALSE]), type="response", 
  	  		  	 				X=d[,-1, drop=FALSE], Y=d[,1])[[1]])
  	   } else{
  	   	 	g1W <- cbind(predict(m, newdata=data.frame(Z=0, A=0, newdata[,-(1:3), drop=FALSE]), type="response"), 
  	  			 		  predict(m, newdata=data.frame(Z=0, A=1, newdata[,-(1:3), drop=FALSE]), type="response"),
  	  		  	 		  predict(m, newdata=data.frame(Z=1, A=0, newdata[,-(1:3), drop=FALSE]), type="response"),	     	     		  predict(m, newdata=data.frame(Z=1, A=1, newdata[,-(1:3), drop=FALSE]), type="response"))
  	   }
  	   colnames(g1W) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")
  	  }
  	}	
} else {
  		type <- "user-supplied values"
  		if(outcome=="Z") {
  			colnames(g1W) <- c("A0", "A1")
  		} else if (outcome=="D"){
  			colnames(g1W) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")[1:ncol(g1W)]
  		}
  	}
  	if(is.null(type)){ type <- class(m)[1]}
  	returnVal <- list(g1W=g1W, coef=coef, type=type)
  	if(type=="SuperLearner"){
			 returnVal$SL.library=SL.library
			 returnVal$coef=m$coef
		}
 return(returnVal)
}


#----------------------- calcParameters -----------------------
# Calculate parameter estimates, CIs and p-values
# arguments:
#   Y - outcome
#   A - binary treatment indicator
#   I.Z - Indicator that Z=z (needed for CDE estimation)
#   Delta - missingness indicator
# 	g1W - values of P(Delta=1|Z,A,W)*P(Z=z|A,W)*P(A=1|W)
# 	g0W - values of P(Delta=1|Z,A,W)*P(Z=z|A,W)*P(A=0|W)
#   Delta =- missingness indicator 
#   Q - nx3 matrix(Q(A,W), Q(1,W), Q(0,W))
#   mu1 - targeted estimate of EY1
#   mu0 - targeted estimate of EY0
#   id - subject id
#   family - "gaussian" or "binomial" (or "poisson") 
# returns:
#   list: (EY1, ATE, RR, OR) 
#	 EY1 - population mean outcome (NULL unless no variation in A)
#	 ATE - additive effect: psi, CI, pvalue
#	 RR - relative risk (NULL if family != binomial): psi, CI, pvalue, log.psi, var.log.psi
#	 OR - odds ratio (NULL if family != binomial): psi, CI, pvalue, log.psi, var.log.psi
#	 IC - list of IC for each parameter (RR and OR on log scale)
#--------------------------------------------------------------
calcParameters <- function(Y,A, I.Z, Delta, g1W, g0W, Q, mu1, mu0, id, family){
	n.id <- length(unique(id))
	Y[is.na(Y)] <- 0
	EY1 <- ATE <- RR <- OR <- NULL
	IC.EY1 <- IC.ATE <- IC.logRR <- IC.logOR <- NULL
	if(length(unique(A))==1){
		pDelta1 <- A*g1W+ (1-A)*g0W   # one of these will be zero, the other pDelta1
		IC.EY1 <- Delta/pDelta1*(Y-Q[,"QAW"]) + Q[,"Q1W"] - mu1
		if(n.id < length(id)){
			IC.EY1 <- as.vector(by(IC.EY1, id, mean))
		}
		EY1$psi <- mu1
		EY1$var.psi <- var(IC.EY1)/n.id
		EY1$CI <- EY1$psi + c(-1.96,1.96)*sqrt(EY1$var.psi)
		EY1$pvalue <- 2*pnorm(-abs(EY1$psi/sqrt(EY1$var.psi)))
	} else {
		IC.ATE<- I.Z*(A/g1W - (1-A)/g0W)*Delta*(Y-Q[,"QAW"]) + Q[,"Q1W"] - Q[,"Q0W"] - (mu1-mu0)
		if(n.id < length(id)){
			IC.ATE <- as.vector(by(IC.ATE, id, mean))
		}
		ATE$psi <- mu1-mu0
		ATE$var.psi <- var(IC.ATE)/n.id
		ATE$CI <- ATE$psi + c(-1.96,1.96)*sqrt(ATE$var.psi)
		ATE$pvalue <- 2*pnorm(-abs(ATE$psi/sqrt(ATE$var.psi)))

		if(family=="binomial"){		
			IC.logRR <- 1/mu1*(I.Z*(A/g1W)*Delta*(Y-Q[,"QAW"]) + Q[,"Q1W"] - mu1) - 
					1/mu0*(I.Z*(1-A)/g0W*Delta*(Y-Q[,"QAW"])+ Q[,"Q0W"] - mu0)
			if(n.id < length(id)){
				IC.logRR <- as.vector(by(IC.logRR, id, mean))
			}
			var.psi.logRR <- var(IC.logRR)/n.id
			RR$psi <- mu1/mu0
			RR$CI  <- exp(log(RR$psi) + c(-1.96,1.96)*sqrt(var.psi.logRR))
			RR$pvalue <- 2*pnorm(-abs(log(RR$psi)/sqrt(var.psi.logRR)))
			RR$log.psi <- log(RR$psi)
			RR$var.log.psi <- var.psi.logRR

			IC.logOR <- 1/(mu1*(1-mu1)) * (I.Z*A/g1W*Delta*(Y-Q[,"QAW"]) + Q[,"Q1W"]) - 
					1/(mu0*(1-mu0)) * (I.Z*(1-A)/g0W*Delta * (Y-Q[,"QAW"]) + Q[,"Q0W"])
			if(n.id < length(id)){
				IC.logOR <- as.vector(by(IC.logOR, id, mean))
			}
			var.psi.logOR <- var(IC.logOR)/n.id
			OR$psi <-  mu1/(1-mu1)/(mu0 / (1-mu0))
			OR$CI  <- exp(log(OR$psi) + c(-1.96,1.96)*sqrt(var.psi.logOR))
			OR$pvalue <- 2*pnorm(-abs(log(OR$psi)/sqrt(var.psi.logOR)))
			OR$log.psi <- log(OR$psi)
			OR$var.log.psi <- var.psi.logOR
		}
	}
	return(list(EY1=EY1, ATE=ATE, RR=RR, OR=OR, IC = list(IC.EY1=IC.EY1, IC.ATE=IC.ATE, IC.logRR = IC.logRR, IC.logOR = IC.logOR)))
}

#-------------------------------tmle----------------------------------------
# estimate marginal treatment effect for binary point treatment
# accounting for missing outcomes.
# EY1 parameter if no variation in A
# arguments:
# Y - outcome
# A - binary treatment indicator, 1-treatment, 0 - control
# W - vector, matrix or dataframe containing baseline covariates
# Z - optional binary intermediate between A and Y - if specified
# 	  we'll compute "controlled direct effect" to get parameter estimate at
#    each value of Z 
# Delta - indicator of missingness.  1 - observed, 0 - missing
# Q - E(Y|Z,A,W), optional nx2 matrix [E(Y|A=0,W), E(Y|A=1,W)] (if CDE estimation, Z=0)
# Q.Z1 - optional nx2 matrix [E(Y|A=0,W), E(Y|A=1,W)] (if CDE estimation, when Z=1)
# g1W - optional values for P(A=1|W)
# gform - optional glm regression formula 
# gbound - one value for symmetrical bounds on g1W, or a vector containing upper and lower bounds
# pZ1 - optional values for P(Z=1|A,W)
# g.Zform - optional glm regression formula
# pDelta1 - optional values for P(Delta=1|Z,A,W)
# g.Deltaform - optional glm regression formula  
# Q.SL.library- optional Super Learner library for estimation of Q, 
# cvQinit - if TRUE obtain cross-validated initial Q
# g.SL.library - optional library for estimation of g and g.Delta
# family - family specification for regression models, defaults to gaussian
# fluctuation - "logistic" (default) or "linear" (for targeting step)
# alpha - bound on predicted probabilities for Q (0.005, 0.995 default)
# id - optional subject identifier
# verbose - flag for controlling printing of messages
#-------------------------------------------------------------------------------
tmle <- function(Y,A,W,Z=NULL, Delta=rep(1,length(Y)),  
				Q=NULL, Q.Z1=NULL, Qform=NULL, Qbounds=NULL, 
				Q.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"), cvQinit=FALSE,
				g1W=NULL, gform=NULL, gbound=0.025, 
				pZ1=NULL, g.Zform=NULL,
				pDelta1=NULL, g.Deltaform=NULL, 
				g.SL.library=c("SL.glm", "SL.step", "SL.glm.interaction"),
				family="gaussian", fluctuation="logistic",
				alpha  = 0.995, id=1:length(Y),verbose=FALSE) {
	# Initializations
	psi.tmle <- varIC <- CI <- pvalue <- NA
	# W <- as.matrix(W)  # if W is a dataframe with factors, changing to matrix will blow it
	colnames(W) <- .setColnames(colnames(W), NCOL(W), "W")
	if (identical(family, binomial)) {
        family == "binomial"
    }else if (identical(family, gaussian)) {
        family == "gaussian"
    }else if (identical(family, poisson)) {  
        family == "poisson"
        if(is.null(Qform)){ 		# only glm for Q when family=poisson
        	Qform <- paste("Y~A", paste(colnames(W), collapse="+"), sep="+")
        }
    }
	if(is.null(A) | all(A==0)){
		A <- rep(1, length(Y))
	}
	
	if(!.verifyArgs(Y,Z,A,W,Delta, Qform, gform, g.Zform, g.Deltaform)){
		stop()
	}
   
 	maptoYstar <- fluctuation=="logistic"
 		    	
   	if(!is.null(Z) & !is.null(pZ1)){
   		if(NCOL(pZ1)==2){
   			colnames(pZ1) <- c("A0", "A1")
   		} else {
   			stop("pZ1 must be an nx2 matrix: [P(Z=1|A=0,W),P(Z=1|A=1,W)]\n")
   		}
   	}
   	
  	if(any(Delta!=1) & !is.null(pDelta1)){
   		if(NCOL(pDelta1)==2 & is.null(Z)){
   			pDelta1 <- cbind(pDelta1, pDelta1)
   			colnames(pDelta1) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")  # won't use the second set
   		} else if(NCOL(pDelta1)==4){
   			colnames(pDelta1) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1") 
   		}else {
   			if(is.null(Z)){
   				stop("pDelta1 must be an nx2 matrix: [P(Delta=1|A=0,W), P(Delta=1|A=1,W)]\n")
   			} else {
   				stop("pDelta1 must be an nx4 matrix:\n  [P(Delta=1|Z=0,A=0,W), P(Delta=1|Z=0,A=1,W), P(Delta=1|Z=1,A=0,W), P(Delta=1|Z=1,A=1,W)]\n")
   			}
   		}
   	}
   	   	
   if(is.null(Z)){ 
   		Z <- rep(1, length(Y))	
   	}
   	CDE <- length(unique(Z))>1
	
	# Stage 1	
 	stage1 <- .initStage1(Y, A, Q, Q.Z1, Delta, Qbounds, alpha, maptoYstar, family)		
	Q <- suppressWarnings(estimateQ(Y=stage1$Ystar,Z,A,W, Delta, Q=stage1$Q, Qbounds=stage1$Qbounds, Qform, 
					maptoYstar=maptoYstar, SL.library=Q.SL.library, 
					cvQinit=cvQinit, family=family, id=id, verbose=verbose))
					
	# Stage 2
	if(length(gbound)==1){
		if(length(unique(A))==1 & length(unique(Z))==1){  # EY1 only, no controlled direct effect
			gbound <- c(gbound,1)
		} else {
			gbound <- c(gbound, 1-gbound)
		}
	}
 	g <- suppressWarnings(estimateG(d=data.frame(A,W), g1W, gform, g.SL.library, id=id, verbose, "treatment mechanism", outcome="A")) 
 	g$bound <- gbound
 	if(g$type=="try-error"){
 		stop("Error estimating treatment mechanism (hint: only numeric variables are allowed)") 
 	}
	
  	if(!CDE){
  		g.z <- NULL
  		g.z$type="No intermediate variable"
  		g.z$coef=NA
  		g.Delta <- suppressWarnings(estimateG(d=data.frame(Delta, Z=1, A, W), pDelta1, g.Deltaform, 
 	 		g.SL.library,id=id, verbose, "missingness mechanism", outcome="D")) 
 		g1W.total <- .bound(g$g1W*g.Delta$g1W[,"Z0A1"], gbound)
  		g0W.total <- .bound((1-g$g1W)*g.Delta$g1W[,"Z0A0"], gbound)  
  		if(all(g1W.total==0)){g1W.total <- rep(10^-9, length(g1W.total))}
  		if(all(g0W.total==0)){g0W.total <- rep(10^-9, length(g0W.total))}
  		H1W <- A/g1W.total
  		H0W <- (1-A)/g0W.total

  		suppressWarnings(
  			epsilon <- coef(glm(stage1$Ystar~-1 + offset(Q$Q[,"QAW"]) + H0W + H1W, family=Q$family, subset=Delta==1))
  		)
  		epsilon[is.na(epsilon)] <- 0  # needed for EY1 calculation
 		Qstar <- Q$Q + c((epsilon[1]*H0W + epsilon[2]*H1W), epsilon[1]/g0W.total, epsilon[2]/g1W.total)
		colnames(Qstar) <- c("QAW", "Q0W", "Q1W")
       Ystar <- stage1$Ystar
   		if (maptoYstar) {
			Qstar <- plogis(Qstar)*diff(stage1$ab)+stage1$ab[1] 
			Q$Q <- plogis(Q$Q)*diff(stage1$ab)+stage1$ab[1]
			Ystar <- Ystar*diff(stage1$ab)+stage1$ab[1]
		} else if (family == "poisson"){  	
    		Q$Q <- exp(Q$Q)				  
    		Qstar <- exp(Qstar)
    	}
    	colnames(Q$Q) <- c("QAW", "Q0W", "Q1W")
    	Q$Q <- Q$Q[,-1]
    	res <- calcParameters(Ystar, A, I.Z=rep(1, length(Ystar)), Delta, g1W.total, g0W.total, Qstar, 
   	   		  		mu1=mean(Qstar[,"Q1W"]), mu0=mean(Qstar[,"Q0W"]), id, family)
  		returnVal <- list(estimates=res, Qinit=Q, g=g, g.Z=g.z, g.Delta=g.Delta, Qstar=Qstar[,-1], epsilon=epsilon) 
  		class(returnVal) <- "tmle"
  	} else {
  		returnVal <- vector(mode="list", length=2)
  		g.z <- suppressWarnings(estimateG(d=data.frame(Z,A,W), pZ1, g.Zform, g.SL.library, id=id, 
  					  verbose, "intermediate variable", outcome="Z"))
  		g.Delta <- suppressWarnings(estimateG(d=data.frame(Delta,Z, A, W), pDelta1, g.Deltaform, 
  								 g.SL.library,id=id, verbose, "missingness mechanism", outcome="D")) 
    	ZAD <- cbind(D1Z0A0 = .bound((1-g$g1W)*(1-g.z$g1W[,"A0"])*g.Delta$g1W[,"Z0A0"], gbound),
  					  D1Z0A1 = .bound(g$g1W*(1-g.z$g1W[,"A1"])*g.Delta$g1W[,"Z0A1"], gbound),
  					  D1Z1A0 = .bound((1-g$g1W)*g.z$g1W[,"A0"]*g.Delta$g1W[,"Z1A0"], gbound),
  					  D1Z1A1 = .bound(g$g1W*g.z$g1W[,"A1"]*g.Delta$g1W[,"Z1A1"], gbound))
  	   adjustZero <- colSums(ZAD)==0
  	   ZAD[,adjustZero] <- 10^-9
  	   	for (z in 0:1){
  	     	H0W <- (1-A)*(Z==z)/ZAD[,z*2+1]
  			H1W <- A*(Z==z) /ZAD[,z*2+2]
  			suppressWarnings(
  				epsilon <- coef(glm(stage1$Ystar~-1 + offset(Q$Q[,"QAW"]) + H0W + H1W, family=Q$family,
  					  subset=(Delta==1 & Z==z)))
  			)  			

  			hCounter <- cbind(1/ZAD[,z*2+1], 1/ZAD[,z*2+2])
 			Qstar <- Q$Q[,c(1, z*2+2, z*2+3)] + c((epsilon[1]*H0W + epsilon[2]*H1W), 
 												        epsilon[1]*hCounter[,1], epsilon[2]*hCounter[,2])
			colnames(Qstar) <- c("QAW", "Q0W", "Q1W") 
			newYstar <- stage1$Ystar     
   			if (maptoYstar) {
				Qstar <- plogis(Qstar)*diff(stage1$ab)+stage1$ab[1] 
				Qinit.return <- plogis(Q$Q[,c(1, z*2+2, z*2+3)])*diff(stage1$ab)+stage1$ab[1]
				newYstar <- stage1$Ystar*diff(stage1$ab)+stage1$ab[1]
			} else if (family == "poisson"){ 
    			Qinit.return <- exp(Q$Q[,c(1, z*2+2, z*2+3)]) 
    			Qstar <- exp(Qstar)
    		}
    		colnames(Qinit.return) <- c("QAW", "Q0W", "Q1W")
    		res <- calcParameters(newYstar, A,I.Z=as.integer(Z==z), Delta, g1W=ZAD[,z*2+2], g0W=ZAD[,z*2+1],  
   	   					  Qstar, mu1=mean(Qstar[,"Q1W"]), mu0=mean(Qstar[,"Q0W"]), id, family)
   	   		Qreturn <- Q
   	   		Qreturn$Q <- Qinit.return[,-1]
   	   		returnVal[[z+1]] <- list(estimates=res, Qinit=Qreturn, g=g, g.Z=g.z, g.Delta=g.Delta, 
   	   									 Qstar=Qstar[,-1], epsilon=epsilon)
  		}
  		class(returnVal[[1]]) <- class(returnVal[[2]]) <- "tmle"
  		class(returnVal) <- "tmle.list"
  	}
  	return(returnVal)
}



#Copyright ©2012. The Regents of the University of California (Regents). All Rights Reserved. Permission to use, copy, modify, and distribute this software and its documentation for #educational, research, and not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the #following two paragraphs appear in all copies, modifications, and distributions. Contact The Office of Technology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA #94720-1620, (510) 643-7201, for commercial licensing opportunities.

#[Created by Susan Gruber, Sam Lendle and Mark van der Laan, Department of Biostatistics, University of California, Berkeley.]

#IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

