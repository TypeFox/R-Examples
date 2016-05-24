#' Load/unload C-lib.

.onLoad <- function(lib, pkg)
{
	library.dynam(chname="VCA", package=pkg, lib.loc=lib)
	# create VCA message environment
	msgEnv <<- new.env(parent=emptyenv())
	check4MKL()
}

.onUnload <- function(lib)
{
	library.dynam.unload(chname="VCA", libpath=lib)
}	


#' Load 'RevoUtilsMath'-package if available.
#' 
#' This function is taken from the Rprofile.site file of Microsoft R Open.
#' It was added to the package namespace to avoid a NOTE during the R CMD check
#' process stating that this function is not gobally defined.
#' 
#' Only change to the original version is a different bracketing scheme to match
#' the one used in the remaining source-code of the package. 
#' 
#' @param package		(character) package name to load, usually this will be package
#' 						'RevoUtilsMath' if available
#' 
#' @author Authors of the Rprofile.site file in Microsoft R Open.

load_if_installed <- function(package) 
{
	if (!identical(system.file(package="RevoUtilsMath"), "")) 
	{
		do.call('library', list(package))
		return(TRUE)
	} 
	else
		return(FALSE)
}



#' Check for Availability of Intel's Math Kernel Library.
#' 
#' Majority of the code is borrowed from the Microsoft R Open Rprofile.site file.
#' In case MKL can be detected this information will be stored in a separate envrionment, which
#' other function know about. If so, an optimized version of function \code{\link{getGB}}
#' will be used which used ordinary matrix-objects instead of matrices defined by the
#' \code{Matrix}-package. This seems to accelerate computation time for large datasets
#' by up to factor 30.
#' 
#' This function is for internal use only and therefore not exported.
#' 
#' @return variable 'MKL' in envir "msgEnv" will be set to TRUE/FALSE
#' 
#' @author 	Authors of the Rprofile.site file in Microsoft R Open,
#' 			Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

check4MKL <- function()
{
	if("msgEnv" %in% ls(.GlobalEnv) && !is.null(msgEnv$MKL))
		return(msgEnv$MKL)
	else
	{
		msgEnv <<- new.env(parent=emptyenv())				
		
		if(Sys.info()["sysname"] == "Darwin")				# Mac OSx
		{		
			assign("MKL", TRUE, envir=msgEnv)				# set MKL to TRUE, although, it may not be installed --> no good method known to check for MKL under MacOS
		} 
		else 												# other operating systems
		{		
			MRO <- FALSE
			try(MRO <- load_if_installed("RevoUtilsMath"), silent=TRUE)			# function only exists in MRO environment
			if(!inherits(MRO, "try-error") && MRO)
				assign("MKL", TRUE, envir=msgEnv)
			else 
				assign("MKL", FALSE, envir=msgEnv)
		}
		
		return(msgEnv$MKL)
	}
}

#' Giesbrecht & Burns Approximation of the Variance-Covariance Matrix of Variance Components.
#'
#' Compute variance covariance matrix of variance components of a linear mixed model
#' via the method stated in Giesbrecht and Burns (1985). 
#' 
#' This function is not intended to be called by users and therefore not exported.
#'
#' @param obj		(object) with list-type structure, e.g. \code{VCA} object fitted by ANOVA
#' 					or a premature \code{VCA} object fitted by REML
#' @param tol		(numeric) values < 'tol' will be considered being equal to zero
#' 
#' @return 	(matrix) corresponding to the Giesbrecht & Burns approximation
#'			of the variance-covariance matrix of variance components
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @seealso \code{\link{vcovVC}}, \code{\link{remlVCA}}, \code{\link{remlMM}}
#'
#' @references
#' Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York
#' 
#' Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#' Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486 

getGB <- function(obj, tol=1e-12)
{
	Nvc	<- obj$Nvc
	Z	<- obj$Matrices$Zre
	Q	<- obj$Matrices$Q
	
	if(is.null(Q))							# should only be missing when fitted by ANOVA
	{
		X   <- getMat(obj, "X")
		Vi  <- getMat(obj, "Vi")
		T   <- getMat(obj, "T")				# P %*% t(X) %*% Vi
		Q   <- Vi - Vi %*% X %*% T			
	}
	VCvar <- matrix(0, Nvc, Nvc)	
	ZxZt  <- vector("list", Nvc)
	nze   <- NULL							# non-zero elements
	
	for(i in 1:Nvc)
	{		
		VCi <- obj$VCoriginal[i]			# for model fitted via ANOVA, decide on basis of orignal (unconstrained) estimates 
		
		if(is.null(VCi))
			VCi <- obj$aov.tab[obj$re.assign$terms[i],"VC"]
		
		if(i < Nvc && abs(VCi) < tol)	
		{
			VCvar[i,] <- VCvar[,i] <- 0 
			next
		}
		else
			nze <- c(nze, i)
		
		for(j in i:Nvc)
		{		
			if(is.null(ZxZt[[i]]))
			{
				if( i == Nvc)				# if all previous i-loops were skipped via "next", this condition-handling is required here for i
				{
					ZxZt[[i]] <- Q 
				}
				else
				{
					Zi  <- Z[, which(obj$re.assign$ind == i)]
					ZiT <- t(Zi) 	
					ZxZt[[i]] <- Zi %*% ZiT %*% Q
				}
			}
			
			if(is.null(ZxZt[[j]]))
			{
				if( j == Nvc)
				{
					ZxZt[[j]] <- Q 
				}
				else
				{
					Zj  <- Z[,which(obj$re.assign$ind == j)]
					ZjT <- t(Zj) 
					ZxZt[[j]] <- Zj %*% ZjT %*% Q
				}
			}
			VCvar[i,j] <- VCvar[j,i] <- round(sum(diag(ZxZt[[i]] %*% ZxZt[[j]])), 12)
		}
	}	
	nze 	<- unique(nze)
	VCnames	<- obj$VCnames
	
	if(!"error" %in% VCnames)					# e.g. missing if 'VarVC=FALSE' in 'remlVCA()'
		VCnames <- c(VCnames, "error")
	
	rownames(VCvar) <- colnames(VCvar) <- VCnames
	VCvar[nze, nze] <- 2 * solve(VCvar[nze, nze])
	VCvar <- VCvar
	
	attr(VCvar, "method") <- "gb"
	VCvar
}


#' Construct Variance-Covariance Matrix of Random Effects for Models Fitted by Function 'lmer'.
#' 
#' This function either restricts the variance-covariance matrix of random effects G to be either
#' diagonal ('cov=FALSE') or to take any non-zero covariances into account (default, 'cov=TRUE').
#' 
#' This function is not intended to be called directly by users and therefore not exported!
#' 
#' @param obj		(object) inheriting from class 'lmerMod'
#' @param cov		(logical) TRUE = in case of non-zero covariances a block diagonal matrix will be constructed,
#'                  FALSE = a diagonal matrix with all off-diagonal element being equal to zero will be contructed
#' 
#' @return (Matrix) representing the variance-covariance structure of random effects G
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

lmerG <- function(obj, cov=FALSE)
{
	stopifnot(inherits(obj, "lmerMod"))
	
	vc  <- VarCorr(obj)
	Nre <- unlist(lapply(lme4::ranef(obj), nrow))					# number of random effects per variance component 
	
	lst <- list()
	for(i in 1:length(Nre))
	{
		tmp.vc <- vc[[i]]
		
		if(!cov && nrow(tmp.vc) > 1)
		{
			tmp.vc <- diag(diag(tmp.vc))
		}
		
		lst <- c(lst, replicate(Nre[i], tmp.vc, simplify=FALSE))	# remove covariances from off-diagonal
		
	}
	G <- bdiag(lst)
	G
}


#' Derive and Compute Matrices for Objects Fitted by Function 'lmer'.
#' 
#' Function derives and computes all matrices required for down-stream
#' analyses of VCA-objects fitted with REML via function \code{\link{lmer}}.
#' 
#' Mixed Model Equations (MME) are solved for fixed and random effects applying the same
#' constraints as in \code{\link{anovaMM}}. 
#' The most elaborate and therefore time consuming part is to prepare all matrices required for 
#' approximating the variance-covariance matrix of variance components (\code{\link{getGB}}.
#' To reduce the computational time, this function tries to optimize object-classes depending
#' on whether Intel's (M)ath (K)ernel (L)ibrary could be loaded or not. MKL appears to be more
#' performant with ordinary matrix-objects, whereas all other computations are perfomred using
#' matrix-representations of the \code{Matrix}-package.
#' 
#' This function is not intended to be called directly by users and therefore not exported.
#' 
#' @param obj		(object) inheriting form 'lmerMod'
#' @param tab		(data.frame) representing the basic VCA-table
#' @param terms		(character) vector used for ordering variance components
#' @param cov		(logical) take non-zero covariances among random effects into account (TRUE) or
#' 					not (FALSE), the latter is the default in this package and also implemented in
#' 					\code{\link{remlVCA}}, \code{\link{anovaVCA}}, and \code{\link{anovaMM}}.
#' @param X			(matrix) design matrix of fixed effects as constructed to meet VCA-package requirements
#' 
#' @return (list), a premature 'VCA' object
#' 
#' @seealso \code{\link{remlVCA}}, \code{\link{remlMM}}
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

lmerMatrices <- function(obj, tab=NULL, terms=NULL, cov=FALSE, X=NULL)
{
	stopifnot(inherits(obj, "lmerMod"))
	
	re.org <- re <- lme4::ranef(obj)	# use lme4's ranef S3-method
	
	VCnam  <- NULL
	reInd  <- list()
	last   <- 0
	count  <- 1
	REnam  <- names(re)
	REnamZ <- NULL
	
	for(i in 1:length(re))				# transform random effects into VCA-compatible structure and derive					
	{									# column-index in random effects design matrix Z for all random effects
		if(ncol(re[[i]]) > 1)			# regression model, i.e. random effects in multi-column matrix
		{
			REi 	<- re[[i]]
			NRi 	<- nrow(REi)
			NCi 	<- ncol(REi)
			trm 	<- names(re)[i]
			tmp.re 	<- ind <- NULL			
			
			for(j in 1:NCi)
			{
				reInd[[count]] <- seq(last+1, last+NRi*NCi, by=NCi)		# columns in Z corresponding to a specific random effect
				
				if(j == 1)
				{
					tmp.nam <- paste(trm, rownames(re[[i]]), sep="")
					names(reInd)[count] <- REnam[count]
				}
				else
				{
					tmp.nam <- c(tmp.nam, paste(colnames(REi)[j], tmp.nam, sep=":"))
					names(reInd)[count] <- paste(colnames(REi)[j], REnam[i], sep=":")
				}
				last  <- last + 1
				count <- count + 1
				
				if(j == 1)
					ind <- eval(parse(text=paste("((",j,"-1)*",NRi,"+1):(", j,"*", NRi, ")", sep="")))
				else
					ind <- paste(ind, eval(parse(text=paste("((",j,"-1)*",NRi,"+1):(", j,"*", NRi, ")", sep=""))), sep=", ")
			}
			ind 	<- paste("c(", paste(ind, collapse=","), ")", sep="")
			ind 	<- eval(parse(text=ind))
			REnamZ 	<- c(REnamZ, tmp.nam[ind])							# column names in Z-matrix			
			cn 		<- colnames(REi)
			cn 		<- paste(names(re)[i], cn, sep=":")
			cn 		<- gsub(":\\(Intercept\\)", "", cn)					# get rif of ()
			VCnam 	<- c(VCnam, cn)										# names of Variance components
		}
		else		# single column random effects matrix
		{
			reInd[[count]] <- (last+1):(last+nrow(re[[i]]))
			names(reInd)[count] <- REnam[count]
			last   <- last + nrow(re[[i]])
			count  <- count + 1
			nami   <- unlist(strsplit(REnam[i], ":"))
			rni    <- rownames(re[[i]])
			VCnam  <- c(VCnam, REnam[i])			
			REnamZ <- c(REnamZ, sapply(rni, function(x) paste(paste(nami, unlist(strsplit(x, ":")), sep=""), collapse=":")))
		}
	}
	REnamZ 		<- as.character(REnamZ)						# names for the random effects and the columns in Z
	reInd  		<- reInd[terms]								# order according to order of variables in user-specified formula
	reInd  		<- reInd[which(!is.na(names(reInd)))]
	Nre	   		<- as.numeric(unlist(lapply(re.org, nrow)))
	sNre   		<- sum(Nre)
	re.assign 	<- list(ind=integer(sNre), terms=names(reInd))
	VC	  	  	<- tab[-c(1, nrow(tab)), "VC"]
	
	for(i in 1:length(reInd))
		re.assign$ind[reInd[[i]]] <- i
	
	Nvc	  <- nrow(tab)-1							# minus total and error								
	Zt 	  <- obj@pp$Zt
	
	if(check4MKL())
		Zt <- as.matrix(Zt)
	
	if(is.null(X) || class(X) != "matrix")
		X <- model.matrix(obj, type="fixed")
	
	Xt <- t(X)
	Z  <- t(Zt)
	colnames(Z) <- REnamZ
	G  <- lmerG(obj, cov=cov)						# construct G-matrix
	
	if(check4MKL())
		R	<- diag(nrow(obj@frame))* tab[nrow(tab),"VC"]
	else
		R	<- Diagonal(nrow(obj@frame))* tab[nrow(tab),"VC"]
	
	V	 <- Z %*% G %*% Zt + R						# variance-covariance matrix of observations

	Vi	 <- solve(V)
	ViX  <- Vi %*% X
	XtVi <- Xt %*% Vi
	P 	 <- MPinv(Xt %*% ViX)   					# re-compute since 'vcov(obj)' differs from e.g. SAS PROC MIXED
	Q	 <- Vi - ViX %*% P %*% XtVi
	T    <- P %*% XtVi
	
	res 		 <- list()
	res$Nvc 	 <- Nvc
	res$VCnames  <- c(terms, "error") 
	y  <- obj@resp$y
	fe <- P %*% XtVi %*% y						# solve mixed model equations for fixed effects
	
	colnames(fe) <- "Estimate"
	rownames(fe) <- colnames(X)
	iind <- which(names(fe) == "(Intercept)")
	
	if(length(iind) > 0)
		names(fe)[iind] <- "int"
	
	res$Matrices 		<- list(Zre=Z, G=G, R=R, V=V, Vi=Vi, Q=Q, X=X, T=T, y=y)
	res$FixedEffects  	<- fe	
	re  				<- G %*% Zt %*% Vi %*% (y - X %*% fe) 		# estimate random effects solving Mixed Model Equations
	re 					<- Matrix(re)								# re-order random effects accordingly
	rownames(re) 		<- REnamZ
	res$RandomEffects 	<- re
	res$re.assign 		<- re.assign
	res$VarFixed  		<- P
	res$ColOrderZ 		<- reInd			# keep information about original col-indexing
	res
}

#' Derive VCA-Summary Table from an object fitted via function \code{\link{lmer}}.
#'
#' This function builds a variance components analysis results table
#' from an object representing a model fitted by \code{\link{lmer}} of the
#' \code{lme4} R-package. It applies the approximation of the variance-covariance
#' matrix of variance components according to Giesbrecht & Burns (1985) and uses this
#' information to approximate the degrees of freedom according to Satterthwaite
#' (see SAS PROC MIXED documentation option 'CL').
#' 
#' This function is not intended to be called directly by users and therefore not exported.
#'
#' @param obj		(lmerMod) object as returned by function lmer
#' @param VarVC		(logical) TRUE = the variance-covariance matrix of variance components will be approximated
#' 					following the Giesbrecht & Burns approach, FALSE = it will not be approximated	
#' @param terms		(character) vector, optionally defining the order of variance terms to be used
#' @param Mean		(numeric) mean value used for CV-calculation
#' @param cov		(logical) TRUE = in case of non-zero covariances a block diagonal matrix will be constructed,
#'                  FALSE = a diagonal matrix with all off-diagonal element being equal to zero will be contructed
#' @param X			(matrix) design matrix of fixed effects as constructed to meet VCA-package requirements
#' 
#' @return (list) still a premature 'VCA' object but close to a
#' 
#' @seealso \code{\link{remlVCA}}, \code{\link{remlMM}}
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#' @references
#' Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York
#' 
#' Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#' Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486 

lmerSummary <- function(obj, VarVC=TRUE, terms=NULL, Mean=NULL, cov=FALSE, X=NULL)
{
	stopifnot(inherits(obj, "lmerMod"))	

	Sum  <- as.data.frame(summary(obj)$varcor)
	Sum[which(Sum[,"var1"] %in% c(NA, "(Intercept)")), "var1"] <- ""
	
	Sum[,"grp"] <- apply(Sum[,c("grp", "var1")], 1, function(x){
				if(x[2] == "")
					return(x[1])
				else
					return(paste(rev(x), collapse=":"))
			})
	re.cor <- NULL
	
	if(any(!is.na(Sum[,"var2"])))
	{
		ind.cor <- which(!is.na(Sum[,"var2"]))
		re.cor  <- Sum[ind.cor,,drop=FALSE]
		Sum <- Sum[-ind.cor,]
	}
	
	Sum  <- Sum[,-c(2,3)]
	rownames(Sum) <- Sum[,"grp"]
	Sum <- Sum[c(terms, "Residual"), ]
	colnames(Sum) <- c("Name", "VC", "SD")
	Sum <- rbind(c(Name="total", VC=sum(Sum[,"VC"]), SD=sqrt(sum(Sum[,"VC"]))), Sum)
	rownames(Sum) <- Sum[,"Name"]
	Sum$VC <- as.numeric(Sum$VC)
	Sum$SD <- as.numeric(Sum$SD)
	Sum$Perc <- 100*Sum$VC/Sum$VC[1]
	Sum$CV	 <- 100*Sum$SD/Mean
	obj <- lmerMatrices(obj, tab=Sum, terms=terms, 		# compute some required matrices
			cov=cov, X=X)	
	obj$aov.tab <- Sum
	
	if(VarVC)
	{
		varVC <- obj$VarCov <- getGB(obj)				# apply Giesbrecht & Burns approximation optimized for MKL
		varVC <- diag(varVC)
		varVC <- c(sum(obj$VarCov), varVC)				# variance of total is sum of all elements of the variance-covariance matrix
		seVC  <- sqrt(varVC)
		Sum$varVC <- varVC
		Sum$Wald <- Sum$VC/seVC							# Wald-statistic
		Sum$DF	 <- 2*Sum$Wald^2						# see SAS-Doc of PROC MIXED
		Sum 	 <- Sum[,c(8,2,4,3,5,6)]
		colnames(Sum)[c(3,5,6)] <- c("%Total", "CV[%]", "Var(VC)")
	}
	else
	{	
		Sum <- Sum[,c(2,4,3,5)]
		colnames(Sum)[c(2,4)] <- c("%Total", "CV[%]")
		obj$aov.tab <- Sum
	}
	
	if(check4MKL())			# remaining part of the package computes with Matrix-package
	{
		obj$Matrices$Zre 	<- Matrix(obj$Matrices$Zre)
		obj$Matrices$G   	<- Matrix(obj$Matrices$G)
		obj$Matrices$R 		<- Matrix(obj$Matrices$R)
		obj$Matrices$V 		<- Matrix(obj$Matrices$V)
		obj$Matrices$Vi 	<- Matrix(obj$Matrices$Vi)
		obj$Matrices$Q 		<- Matrix(obj$Matrices$Q)
		obj$Matrices$X	    <- Matrix(obj$Matrices$X)
		obj$Matrices$T	    <- Matrix(obj$Matrices$T)
	}	
	
	rownames(Sum)[nrow(Sum)] <- "error"
	obj$aov.tab <- Sum
	obj$re.cor <- re.cor					# save correlation among random terms
	obj
}


#' Perform (V)ariance (C)omponent (A)nalysis via REML-Estimation.
#' 
#' Function performs a Variance Component Analysis (VCA) using Restricted Maximum Likelihood (REML)
#' to fit the random model, i.e. a linear mixed model (LMM) where the intercept is the only fixed effect.
#' 
#' Here, a variance component model is fitted by REML using the \code{\link{lmer}} function of the
#' \code{lme4}-package. For all models the Giesbrechnt & Burns (1985) approximation of the variance-covariance
#' matrix of variance components (VC) is applied. A Satterthwaite approximation of the degrees of freedom
#' for all VC and total variance is based on this approximated matrix using \eqn{df=2Z^2}, where
#' \eqn{Z} is the Wald statistic \eqn{Z=\sigma^2/se(\sigma^2)}, and \eqn{\sigma^2} is here used for an
#' estimated variance. The variance of total variability, i.e. the sum of all VC is computed via summing
#' up all elements of the variance-covariance matrix of the VC.
#' Note, that for large datasets approximating the variance-covariance matrix of VC is computationally expensive
#' and may take very long. There is no Fisher-information matrix available for 'merMod' objects, which can
#' serve as approximation. To avoid this time-consuming step, use argument 'VarVC=FALSE' but remember,
#' that no confidence intervals for any VC will be available. If you use Microsoft's R Open, formerly known
#' as Revolution-R, which comes with Intel's Math Kernel Library (MKL), this will be automatically detected
#' and an environment-optimized version will be used, reducing the computational time very much (see examples).
#' 
#' @param form          (formula) specifying the model to be fit, a response variable left of the '~' is mandatory
#' @param Data          (data.frame) storing all variables referenced in 'form'
#' @param by			(factor, character) variable specifying groups for which the analysis should be performed individually,
#' 						i.e. by-processing
#' @param VarVC			(logical) TRUE = the variance-covariance matrix of variance components will be approximated using 
#' 						the method found in Giesbrecht & Burns (1985), which also serves as basis for applying a Satterthwaite
#' 						approximation of the degrees of freedom for each variance component, FALSE = leaves out this step, 
#' 						no confidence intervals for VC will be available
#' @param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' 
#' @seealso \code{\link{remlMM}}, \code{\link{VCAinference}}, \code{\link{ranef.VCA}}, \code{\link{residuals.VCA}},
#' 			\code{\link{anovaVCA}}, \code{\link{anovaMM}}, \code{\link{plotRandVar}}, \code{\link{lmer}}
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples
#' \dontrun{
#' 
#' # a VCA standard example
#' data(dataEP05A2_3)
#' 
#' # fit it by ANOVA first, then by REML
#' fit0 <- anovaVCA(y~day/run, dataEP05A2_3) 
#' fit1 <- remlVCA(y~day/run, dataEP05A2_3)
#' fit0
#' fit1
#'  
#' # make example unbalanced
#' set.seed(107)
#' dat.ub <- dataEP05A2_3[-sample(1:80, 7),]
#' fit0ub <- anovaVCA(y~day/run, dat.ub) 
#' fit1ub <- remlVCA(y~day/run, dat.ub) 
#' 
#' # not that ANOVA- and REML-results now differ
#' fit0ub
#' fit1ub
#' 
#' ### Use the six sample reproducibility data from CLSI EP5-A3
#' ### and fit per sample reproducibility model
#' data(CA19_9)
#' fit.all <- remlVCA(result~site/day, CA19_9, by="sample")
#' 
#' reproMat <- data.frame(
#'  Sample=c("P1", "P2", "Q3", "Q4", "P5", "Q6"),
#'  Mean= c(fit.all[[1]]$Mean, fit.all[[2]]$Mean, fit.all[[3]]$Mean, 
#' 	        fit.all[[4]]$Mean, fit.all[[5]]$Mean, fit.all[[6]]$Mean),
#' 	Rep_SD=c(fit.all[[1]]$aov.tab["error","SD"], fit.all[[2]]$aov.tab["error","SD"],
#' 	         fit.all[[3]]$aov.tab["error","SD"], fit.all[[4]]$aov.tab["error","SD"],
#'           fit.all[[5]]$aov.tab["error","SD"], fit.all[[6]]$aov.tab["error","SD"]),
#' 	Rep_CV=c(fit.all[[1]]$aov.tab["error","CV[%]"],fit.all[[2]]$aov.tab["error","CV[%]"],
#'           fit.all[[3]]$aov.tab["error","CV[%]"],fit.all[[4]]$aov.tab["error","CV[%]"],
#'           fit.all[[5]]$aov.tab["error","CV[%]"],fit.all[[6]]$aov.tab["error","CV[%]"]),
#'  WLP_SD=c(sqrt(sum(fit.all[[1]]$aov.tab[3:4,"VC"])),sqrt(sum(fit.all[[2]]$aov.tab[3:4, "VC"])),
#'           sqrt(sum(fit.all[[3]]$aov.tab[3:4,"VC"])),sqrt(sum(fit.all[[4]]$aov.tab[3:4, "VC"])),
#'           sqrt(sum(fit.all[[5]]$aov.tab[3:4,"VC"])),sqrt(sum(fit.all[[6]]$aov.tab[3:4, "VC"]))),
#'  WLP_CV=c(sqrt(sum(fit.all[[1]]$aov.tab[3:4,"VC"]))/fit.all[[1]]$Mean*100,
#'           sqrt(sum(fit.all[[2]]$aov.tab[3:4,"VC"]))/fit.all[[2]]$Mean*100,
#'           sqrt(sum(fit.all[[3]]$aov.tab[3:4,"VC"]))/fit.all[[3]]$Mean*100,
#'           sqrt(sum(fit.all[[4]]$aov.tab[3:4,"VC"]))/fit.all[[4]]$Mean*100,
#'           sqrt(sum(fit.all[[5]]$aov.tab[3:4,"VC"]))/fit.all[[5]]$Mean*100,
#'           sqrt(sum(fit.all[[6]]$aov.tab[3:4,"VC"]))/fit.all[[6]]$Mean*100),
#'  Repro_SD=c(fit.all[[1]]$aov.tab["total","SD"],fit.all[[2]]$aov.tab["total","SD"],
#'             fit.all[[3]]$aov.tab["total","SD"],fit.all[[4]]$aov.tab["total","SD"],
#'             fit.all[[5]]$aov.tab["total","SD"],fit.all[[6]]$aov.tab["total","SD"]),
#'  Repro_CV=c(fit.all[[1]]$aov.tab["total","CV[%]"],fit.all[[2]]$aov.tab["total","CV[%]"],
#'             fit.all[[3]]$aov.tab["total","CV[%]"],fit.all[[4]]$aov.tab["total","CV[%]"],
#'             fit.all[[5]]$aov.tab["total","CV[%]"],fit.all[[6]]$aov.tab["total","CV[%]"]))
#'  
#'  for(i in 3:8) reproMat[,i] <- round(reproMat[,i],digits=ifelse(i%%2==0,1,3))
#'  reproMat
#' 
#' # now plot the precision profile over all samples
#' plot(reproMat[,"Mean"], reproMat[,"Rep_CV"], type="l", main="Precision Profile CA19-9",
#' 		xlab="Mean CA19-9 Value", ylab="CV[%]")
#' grid()
#' points(reproMat[,"Mean"], reproMat[,"Rep_CV"], pch=16)
#' 
#' 
#' # REML-estimation not yes optimzed to the same degree as
#' # ANOVA-estimation. Note, that no variance-covariance matrix
#' # for the REML-fit is computed (VarVC=FALSE)!
#' # Note: A correct analysis would be done per-sample, this is just
#' #       for illustration.
#' data(VCAdata1)
#' system.time(fit0 <- anovaVCA(y~sample+(device+lot)/day/run, VCAdata1))
#' system.time(fit1 <- remlVCA(y~sample+(device+lot)/day/run, VCAdata1, VarVC=FALSE))
#' 
#' # The previous example will also be interesting for environments using MKL.
#' # Run it once in a GNU-R environment and once in a MKL-environment
#' # and compare computational time of both. Note, that 'VarVC' is now set to TRUE
#' # and variable "sample" is put into the brackets increasing the number of random
#' # effects by factor 10. On my Intel Xeon E5-2687W 3.1 GHz workstation it takes
#' # ~ 400s with GNU-R and ~25s with MKL support (MRO) both run under Windows.
#' system.time(fit2 <- remlVCA(y~(sample+device+lot)/day/run, VCAdata1, VarVC=TRUE))
#' 
#' # using the SWEEP-Operator is even faster but the variance-covariance matrix of
#' # VC is not automatically approximated as for fitting via REML
#' system.time(fit3 <- anovaVCA(y~(sample+device+lot)/day/run, VCAdata1))
#' fit2
#' fit3
#' }

remlVCA <- function(form, Data, by=NULL, VarVC=TRUE, quiet=FALSE)
{
	Call <- match.call()
	
	if(!is.null(by))
	{
		stopifnot(is.character(by))
		stopifnot(by %in% colnames(Data))
		stopifnot(is.factor(by) || is.character(by))
		
		levels  <- unique(Data[,by])
		res <- lapply(levels, function(x) remlVCA(form=form, Data[Data[,by] == x,], VarVC=VarVC, quiet=quiet))
		names(res) <- paste(by, levels, sep=".")
		return(res)
	}
	
	stopifnot(class(form) == "formula")
	stopifnot(is.logical(VarVC))
	stopifnot(is.logical(quiet))
	stopifnot(is.data.frame(Data))
	stopifnot(nrow(Data) > 2)                                               # at least 2 observations for estimating a variance

	trms <- terms(form)														# convert VCA-formula to valid lmer-formula
	stopifnot(attr(trms, "response") == 1)
	lab  <- attr(trms, "term.labels")
	if(length(lab) == 0)
	{
		if(!quiet)
			warning("No random effects specified! Call function 'anovaVCA' instead!")
		return(anovaVCA(form, Data))
	}

	lab  <- paste("(1|", lab, ")", sep="")
	resp <- rownames(attr(trms, "factors"))[1]
	vars <- rownames(attr(trms, "factors"))[-1]
		
	form <- as.formula(paste(resp, "~", paste(lab, collapse="+"), sep=""))
	
	stopifnot(resp %in% colnames(Data))
	stopifnot(is.numeric(Data[,resp]))
	
	rmInd 	<- integer()	
	resp.NA <- is.na(Data[,resp])
	
	if(any(resp.NA))
	{    
		rmInd <- c(rmInd, which(resp.NA))
		if(!quiet)
			warning("There are ", length(which(resp.NA))," missing values for the response variable (obs: ", paste(which(resp.NA), collapse=", "), ")!")
	}    
	
	if(!is.null(vars))
	{
		for(i in vars)                                                  # convert all nested factors as factor-objects
		{
			if( any(is.na(Data[,i])))
			{
				NAind <- which(is.na(Data[,i]))
				rmInd <- c(rmInd, NAind)
				if(!quiet)
					warning("Variable '", i,"' has ",length(NAind)," missing values (obs: ", paste(NAind, collapse=", "), ")!" )
			}
			
			if(length(levels(Data[,i])) >= 0.75*nrow(Data) && !quiet)
				warning("Variable >>> ", i," <<< has at least 0.75 * nrow(Data) levels!")
		}
		rmInd <- unique(rmInd)
	}	
	
	for(i in rev(vars))													# order data
	{
		tmp <- Data[,i]
		tmp.num <- suppressWarnings(as.numeric(as.character(tmp)))		# warnings are very likely here
		
		if(!any(is.na(tmp.num)))										# all levels are numbers
		{
			Width <- max(nchar(as.character(unique(Data[,i])))) + 1
			tmp   <- suppressWarnings(formatC(Data[,i], width=Width))
			uLvl  <- unique(Data[,i])
			Data  <- Data[order(tmp),]
			Data[,i] <- factor(as.character(Data[,i]), levels=uLvl[order(unique(tmp))])
		}		
		else
			Data[,i] <- factor(Data[,i])								# automatically orders
	}
	
	vcol <- rownames(attr(trms, "factors"))
	if(is.null(vcol))
		vcol <- resp
	Ndata <- nrow(Data)
	Data  <- na.omit(Data[,vcol, drop=F])
	Nobs <- nrow(Data)
	
	fit <- lmer(form, Data)													# fit via 'lmer'
	
	res <- list(call=Call, Type="Random Model", EstMethod="REML", data=Data, terms=trms,
			intercept=as.logical(attr(trms, "intercept")), response=resp)
	
	res$Mean 	 	 <- mean(Data[,resp], na.rm=TRUE)
	res$form 	 	 <- form
	res$Nvc  	 	 <- length(lab) + 1											# error is one additional VC
	res$VCnames  	 <- c(attr(trms, "term.labels", "error"))
	res$NegVCmsg 	 <- ""														# there will never be anything to report
	res$VarVC.method <- "gb"
	res$balanced <- if(isBalanced(as.formula(trms), Data)) 
				"balanced"  
			else 
				"unbalanced" 
	
	if(Nobs != Ndata)
		res$Nrm <- Ndata - Nobs                         # save number of observations that were removed due to missing data
	
	res$Nobs <- Nobs
	
	tmp <- lmerSummary(	obj=fit, VarVC=VarVC, 			# construct table similar to aov-table and approximate vcovVC
			terms=attr(trms, "term.labels"),
			Mean=res$Mean)	
	res <- c(res, tmp)
	class(res) <- "VCA"
	res
}


#' Fit Linear Mixed Models via REML.
#' 
#' Function fits Linear Mixed Models (LMM) using Restricted Maximum Likelihood (REML).
#' 
#' Here, a LMM is fitted by REML using the \code{\link{lmer}} function of the \code{lme4}-package. 
#' For all models the Giesbrechnt & Burns (1985) approximation of the variance-covariance
#' matrix of variance components (VC) can be applied ('VarVC=TRUE'). A Satterthwaite approximation of the degrees of freedom
#' for all VC and total variance is based on this approximated matrix using \eqn{df=2Z^2}, where
#' \eqn{Z} is the Wald statistic \eqn{Z=\sigma^2/se(\sigma^2)}, and \eqn{\sigma^2} is here used for an
#' estimated variance. The variance of total variability, i.e. the sum of all VC is computed via summing
#' up all elements of the variance-covariance matrix of the VC.
#' One can constrain the variance-covariance matrix of random effects \eqn{G} to be either diagonal ('cov=FALSE'), i.e.
#' all random effects are indpendent of each other (covariance is 0). If 'cov=TRUE' (the default) matrix \eqn{G} will be
#' constructed as implied by the model returned by function \code{\link{lmer}}. 
#' 
#' As for objects returned by function \code{\link{anovaMM}} linear hypotheses of fixed effects or LS Means can be
#' tested with functions \code{\link{test.fixef}} and \code{\link{test.lsmeans}}. Note, that option "contain" does
#' not work for LMM fitted via REML.
#' 
#' Note, that for large datasets approximating the variance-covariance matrix of VC is computationally expensive
#' and may take very long. There is no Fisher-information matrix available for 'merMod' objects, which can
#' serve as approximation. To avoid this time-consuming step, use argument 'VarVC=FALSE' but remember,
#' that no confidence intervals for any VC will be available. If you use Microsoft's R Open, formerly known
#' as Revolution-R, which comes with Intel's Math Kernel Library (MKL), this will be automatically detected
#' and an environment-optimized version will be used, reducing the computational time considerably (see examples).
#' 
#' @param form          (formula) specifying the model to be fit, a response variable left of the '~' is mandatory
#' @param Data          (data.frame) storing all variables referenced in 'form'
#' @param by			(factor, character) variable specifying groups for which the analysis should be performed individually,
#' 						i.e. by-processing
#' @param VarVC			(logical) TRUE = the variance-covariance matrix of variance components will be approximated using 
#' 						the method found in Giesbrecht & Burns (1985), which also serves as basis for applying a Satterthwaite
#' 						approximation of the degrees of freedom for each variance component, FALSE = leaves out this step, 
#' 						no confidence intervals for VC will be available
#' @param cov			(logical) TRUE = in case of non-zero covariances a block diagonal matrix will be constructed,
#'                  	FALSE = a diagonal matrix with all off-diagonal element being equal to zero will be contructed
#' @param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' 
#' @seealso \code{\link{remlVCA}}, \code{\link{VCAinference}}, \code{\link{ranef.VCA}}, \code{\link{residuals.VCA}},
#' 			\code{\link{anovaVCA}}, \code{\link{anovaMM}}, \code{\link{plotRandVar}},  \code{\link{test.fixef}},  
#' 			\code{\link{test.lsmeans}}, \code{\link{lmer}}
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples
#' \dontrun{
#' data(dataEP05A2_2)
#' 
#' # assuming 'day' as fixed, 'run' as random
#' remlMM(y~day/(run), dataEP05A2_2)
#' 
#' # assuming both as random leads to same results as
#' # calling anovaVCA
#' remlMM(y~(day)/(run), dataEP05A2_2)
#' anovaVCA(y~day/run, dataEP05A2_2)
#' 
#' # use different approaches to estimating the covariance of 
#' # variance components (covariance parameters)
#' dat.ub <- dataEP05A2_2[-c(11,12,23,32,40,41,42),]			# get unbalanced data
#' m1.ub <- remlMM(y~day/(run), dat.ub, SSQ.method="qf", VarVC.method="scm")
#' m2.ub <- remlMM(y~day/(run), dat.ub, SSQ.method="qf", VarVC.method="gb")		# is faster
#' V1.ub <- round(vcovVC(m1.ub), 12)
#' V2.ub <- round(vcovVC(m2.ub), 12)
#' all(V1.ub == V2.ub)
#' 
#' # make it explicit that "gb" is faster than "scm"
#' # compute variance-covariance matrix of VCs 10-times
#' 
#' system.time(for(i in 1:500) vcovVC(m1.ub))	# "scm"
#' system.time(for(i in 1:500) vcovVC(m2.ub))	# "gb"
#' 
#' 
#' # fit a larger random model
#' data(VCAdata1)
#' fitMM1 <- remlMM(y~((lot)+(device))/(day)/(run), VCAdata1[VCAdata1$sample==1,])
#' fitMM1
#' # now use function tailored for random models
#' fitRM1 <- anovaVCA(y~(lot+device)/day/run, VCAdata1[VCAdata1$sample==1,])
#' fitRM1
#' 
#' # there are only 3 lots, take 'lot' as fixed 
#' fitMM2 <- remlMM(y~(lot+(device))/(day)/(run), VCAdata1[VCAdata1$sample==2,])
#' 
#' # the following model definition is equivalent to the one above,
#' # since a single random term in an interaction makes the interaction
#' # random (see the 3rd reference for details on this topic)
#' fitMM3 <- remlMM(y~(lot+(device))/day/run, VCAdata1[VCAdata1$sample==2,])
#' 
#' # fit same model for each sample using by-processing
#' lst <- remlMM(y~(lot+(device))/day/run, VCAdata1, by="sample")
#' lst
#' 
#' # fit mixed model originally from 'nlme' package
#'  
#' library(nlme)
#' data(Orthodont)
#' fit.lme <- lme(distance~Sex*I(age-11), random=~I(age-11)|Subject, Orthodont) 
#' 
#' # re-organize data for using 'remlMM'
#' Ortho <- Orthodont
#' Ortho$age2 <- Ortho$age - 11
#' Ortho$Subject <- factor(as.character(Ortho$Subject))
#' fit.remlMM1 <- remlMM(distance~Sex*age2+(Subject)*age2, Ortho)
#' 
#' # use simplified formula avoiding unnecessary terms
#' fit.remlMM2 <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
#' 
#' # and exclude intercept
#' fit.remlMM3 <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho)
#' 
#' # now use exclude covariance of per-subject intercept and slope
#' # as for models fitted by function 'anovaMM'
#' fit.remlMM4 <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho, cov=FALSE)
#' 
#' # compare results
#' fit.lme
#' fit.remlMM1
#' fit.remlMM2
#' fit.remlMM3
#' fit.remlMM4
#' 
#' # are there a sex-specific differences?
#' cmat <- getL(fit.remlMM3, c("SexMale-SexFemale", "SexMale:age2-SexFemale:age2")) 
#' cmat
#' 			 
#' test.fixef(fit.remlMM3, L=cmat)
#' }

remlMM <- function(form, Data, by=NULL, VarVC=TRUE, cov=TRUE, quiet=FALSE)
{
	Call <- match.call()
	
	if(!is.null(by))
	{
		stopifnot(is.character(by))
		stopifnot(by %in% colnames(Data))
		stopifnot(is.factor(by) || is.character(by))
		
		levels  <- unique(Data[,by])
		res <- lapply(levels, function(x) remlMM(form=form, Data[Data[,by] == x,], VarVC=VarVC, cov=cov, quiet=quiet))
		names(res) <- paste(by, levels, sep=".")
		return(res)
	}
	
	stopifnot(class(form) == "formula")
	stopifnot(is.logical(VarVC))
	stopifnot(is.logical(quiet))
	stopifnot(is.data.frame(Data))
	stopifnot(nrow(Data) > 2)                                               # at least 2 observations for estimating a variance
	
	trms <- terms(form, simplify=TRUE, keep.order=TRUE)						# convert VCA-formula to valid lmer-formula
	stopifnot(attr(trms, "response") == 1)
	resp <- rownames(attr(trms, "factors"))[1]
	org.form <- form
	
	stopifnot(resp %in% colnames(Data))
	stopifnot(is.numeric(Data[,resp]))
	
	res <- list(call=Call,  EstMethod="REML", data=Data, terms=trms,
			response=resp)
	
	int <- res$intercept <- attr(trms, "intercept") == 1						# has intercept	
	rf  <- gregexpr("\\([[:alnum:]]*\\)", as.character(org.form)[3])			# check for random effects
	
	if(rf[[1]][1] != -1)														# identify random variables
	{
		len <- attr(rf[[1]], "match.length")
		pos <- rf[[1]]
		tmp <- NULL
		for(i in 1:length(len))
		{
			tmp <- c(tmp, substr(as.character(org.form)[3], pos[i]+1, pos[i]+len[i]-2))				# remember random factors, exclude brackets
		}
		rf <- tmp
	}
	
	Ndata <- nrow(Data)	
	rmInd <- integer()	
	resp.NA <- is.na(Data[,resp])
	
	if(any(resp.NA))								# remove missing data and warn
	{    
		rmInd <- c(rmInd, which(resp.NA))
		if(!quiet)
			warning("There are ", length(which(resp.NA))," missing values for the response variable (obs: ", paste(which(resp.NA), collapse=", "), ")!")
		res$resp.NA <- rmInd
	}    
	
	fac  	<- attr(trms, "term.labels")
	if(length(fac) > 1)
	{
		rf.ind  <- which(apply(sapply(rf, function(x) regexpr(x, fac)), 1, function(x) any(x>0)))		
	}
	else
	{
		if(length(rf) > 0)
		{
			if(rf == fac)
				rf.ind <- 1
			else
				rf.ind <- numeric(0)
		}
	}
	vars    <- rownames(attr(trms, "factors"))[-1]	                        # remove response
	Nvc     <- length(fac) + 1 
	Order   <- NULL
	
	for(i in rev(vars))														# check Data for consistency
	{
		if( any(is.na(Data[,i])))
		{
			NAind <- which(is.na(Data[,i]))
			rmInd <- c(rmInd, NAind)
			if(!quiet)
				warning("Variable '", i,"' has ",length(NAind)," missing values (obs: ", paste(NAind, collapse=", "), ")!" )
		}
		if(!class(Data[,i]) %in% c("numeric", "character", "factor"))
			stop("Variable ",i," is neither one of \"numeric\", \"character\" or \"factor\"!")
		
		if(class(Data[,i]) != "numeric")
		{
			if(!quiet && is.character(Data[,i]))
				warning("Convert variable ", i," from \"charater\" to \"factor\"!")
			
			tmp <- Data[,i]
			tmp.num <- suppressWarnings(as.numeric(as.character(tmp)))		# warnings are very likely here
			
			if(!any(is.na(tmp.num)))										# all levels are numbers
			{
				Width <- max(nchar(as.character(unique(Data[,i])))) + 1
				tmp   <- suppressWarnings(formatC(Data[,i], width=Width))
				uLvl  <- unique(Data[,i])
				Data  <- Data[order(tmp),]
				Data[,i] <- factor(as.character(Data[,i]), levels=uLvl[order(unique(tmp))])
			}		
			else
				Data[,i] <- factor(Data[,i])								# automatically orders
		}
	}

	if(length(rf.ind) > 0)													# at least one random term in 'form'
	{
		res$random <- fac[rf.ind]
		res$fixed  <- fac[-rf.ind]
	}
	else
	{
		if(!quiet)
			warning("No random terms in the model! Call 'anovaMM' insead!")
		
		return(anovaMM(form, Data))
	}

	res$terms.classes <- sapply(Data[,rownames(attr(trms, "factors"))[-1]], class)
	
	res$Type <- if(length(res$fixed) == 0)
				"Random Model"
			else
				"Mixed Model"			
	
	lab       <- attr(trms, "term.labels")	
	fixed     <- paste(res$fixed, collapse="+")	
	var.num   <- names(res$terms.classes[which(res$terms.classes == "numeric")])
	fac 	  <- attr(trms, "factors")[var.num,res$fixed,drop=FALSE]		# restrict to columns of fixed effects variables and rows with numeric variables
	fe.num 	  <- character()

	fac 	  <- fac[which(apply(fac, 1, any)),,drop=FALSE]
	fac 	  <- fac[,which(apply(fac, 2, any)),drop=FALSE]
	fe.num 	  <- colnames(fac)
	iacs      <- which(grepl(":", res$random))										# interaction terms in the formula
	num.rand  <- sapply(var.num, function(x) grepl(x, res$random)) 
	num.rand  <- as.matrix(num.rand)
	num.fixed <- sapply(var.num, function(x) grepl(x, res$fixed))
	random 	  <- ""
	trmMat    <- matrix(nrow=length(res$random), ncol=2, dimnames=list(NULL, c("form", "subject")))

	for(i in 1:length(res$random))			# process each random term
	{
		if(nchar(random) == 0)
			sep <- ""
		else
			sep <- "+"
		
		if(length(num.rand) > 0 && any(num.rand[i,]))				# random term with numeric variable found
		{
			if(i %in% iacs)
			{
				splt 	<- unlist(strsplit(res$random[i], ":"))
				tmp.num <- which(splt %in% var.num)
				numVar 	<- splt[ tmp.num]
				others  <- paste(splt[-tmp.num], collapse=":")				
				ind 	<- which(trmMat[,2] == others)				# subject terms may only occur once
				
				if(length(ind) > 0)
					trmMat[ind,2] <- NA
				
				trmMat[i,] <- c(numVar, others)	
			}
			else
				stop("Numeric variables may only occur in random interaction terms!")
		}
		else
			trmMat[i,] <- c("1", res$random[i])
	}
	trmMat 	<- na.omit(trmMat)
	random 	<- paste( apply(trmMat, 1, function(x) paste("(", x[1], "|", x[2], ")", sep="")), collapse="+")	
	form 	<- paste(  resp, "~", fixed, "+", random, sep="")
	form 	<- as.formula(form)
	vcol 	<- rownames(attr(trms, "factors"))
	
	if(is.null(vcol))
		vcol <- resp
	
	Ndata <- nrow(Data)
	Data  <- na.omit(Data[,vcol, drop=F])
	Nobs  <- nrow(Data)
	
	fit <- lmer(form, Data)								# fit via 'lmer'
	
	res$Mean 		 <- mean(Data[,resp], na.rm=TRUE)
	res$form 	 	 <- form
	res$NegVCmsg 	 <- ""
	res$VarVC.method <- "gb"
	
	if(int)
	{
		X <- matrix(1, ncol=1, nrow=nrow(Data))			# design matrix of fixed effects: include intercept --> needs a restriction
		colnames(X) <- "int"
		fe.assign <- 0									# '0' indicates intercept
	}
	else
	{
		fe.assign <- NULL
		X <- matrix(1, ncol=0, nrow=nrow(Data))	
	}
	
	fixed.form <- as.formula(paste(	resp, "~", 
									paste(	ifelse(int, "1", "-1"), 
											fixed, sep=ifelse(fixed=="", "", "+")),
									sep=""))
		
	old.opt <- options(contrasts=c("contr.SAS", "contr.poly"))

	suppressWarnings({
		X1 	<- model.matrix(fixed.form, data=Data)						# with SAS contrasts but without the columns where restrictions apply
		X10	<- apply(X1, 2, function(x) all(x==0))
		X2 <- model.matrix(	fixed.form, data=Data,						# full model matrix with re-setted contrasts
							contrasts.arg=lapply(Data[,sapply(Data, is.factor), drop=FALSE],
							contrasts, contrasts=FALSE))
		X20	<- apply(X2, 2, function(x) all(x==0))
		
		X2.asgn 	<- attr(X2, "assign")									# keep info		
		
		if(any(X10))
			X1 <- X1[,-which(X10),drop=FALSE]
		
		if(any(X20))
		{
			X2 <- X2[,-which(X20),drop=FALSE]
			X2.asgn 	<- X2.asgn[-which(X20)]
		}				
	})
	
	options(old.opt)													# reset contrasts option
	
	fixed.terms <- terms(fixed.form)
	fe.terms 	<- attr(fixed.terms, "term.labels")
	
	X2[,!colnames(X2) %in% colnames(X1)] <- 0							# transfer restriction from X1 to X2
	if(int)
	{
		colnames(X2)[1] <- "int"
		fe.terms <- c("int", fe.terms)
	}
	
	fe.assign <- X2.asgn
	attr(fe.assign, "terms") <- fe.terms

	res$fe.assign 	<- fe.assign										# mapping columns of X to fixed terms in the model formula
	res$fixed.terms <- fixed.terms
	res$balanced <- if(isBalanced(as.formula(trms), Data)) 
				"balanced"  
			else 
				"unbalanced" 
	
	if(Nobs != Ndata)
		res$Nrm <- Ndata - Nobs                         # save number of observations that were removed due to missing data
	
	res$Nobs <- Nobs
	
	tmp <- lmerSummary(	obj=fit, VarVC=VarVC, 			# construct table similar to aov-table and approximate vcovVC
						terms=res$random,
						Mean=res$Mean,
						cov=cov, X=X2)				
	
	tmp$Matrices$y <- Matrix(Data[,resp], ncol=1)
	res <- c(res, tmp)
	class(res) <- "VCA"
	res
}


#' Solve System of Linear Equations using Inverse of Cholesky-Root.
#' 
#' This function is intended to reduce the computational time in function
#' \code{\link{solveMME}} which computes the inverse of the square variance-
#' covariance Matrix of observations. It is considerably faster than function
#' \code{\link{solve}} (see example).
#' Whenever an error occurs, which is the case for non positive definite matrices
#' 'X', function \code{\link{MPinv}} is called automatically yielding a generalized
#' inverse (Moore-Penrose inverse) of 'X'.
#' 
#' @param X			(matrix, Matrix) object to be inverted
#' @param quiet		(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' 
#' @return (matrix, Matrix) corresponding to the inverse of X
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples
#' \dontrun{
#' # following complex (nonsense) model takes pretty long to fit
#' system.time(res.sw <- anovaVCA(y~(sample+lot+device)/day/run, VCAdata1))
#' # solve mixed model equations (not automatically done to be more efficient)
#' system.time(res.sw <- solveMME(res.sw))
#' # extract covariance matrix of observations V
#' V1 <- getMat(res.sw, "V")
#' V2 <- as.matrix(V1)
#' system.time(V2i <- solve(V2))
#' system.time(V1i <- VCA:::Solve(V1))
#' V1i <- as.matrix(V1i)
#' dimnames(V1i) <- NULL
#' dimnames(V2i) <- NULL
#' all.equal(V1i, V2i)
#' } 

Solve <- function(X, quiet=FALSE)
{
	stopifnot(ncol(X) == nrow(X))
	
	clsMatrix <- inherits(X, "Matrix")
	if(clsMatrix)
	{
		cls <- class(X)
		X <- as.matrix(X)
	}
	Xi <- try(chol2inv(chol(X)), silent=TRUE)
	
	if(class(Xi) == "try-error")			# use Moore-Penrose inverse instead in case of an error
	{										# using the Cholesky-decomposition approach
		if(!quiet)
			warning("Error in 'chol2inv'!\n\tUse generalized (Moore-Penrose) inverse (MPinv)!", sep="\n")
		Xi <- MPinv(X)
	}
	
	if(clsMatrix)
	{
		Xi <- as(Xi, "dgCMatrix")			# Solve used here for inverting V-matrix
	}
	return(Xi)
}



#' Calling C-implementation of the SWEEP-Operator
#' 
#' Function calls a fast C-implementation of the SWEEP operator using the
#' transpose of the original augmented matrix \eqn{X'X} (see \code{\link{getSSQsweep}}).
#' Transposing prior to applying the SWEEP-operator speeds up things since the 
#' complete matrix is stored in memory in consecutive manner. 
#' 
#' This is an utility-function not intended to be called directly.
#' 
#' @param M			(matrix) matrix, representing the augmented matrix \eqn{X'X}
#' @param asgn		(integer) vector, identifying columns in \eqn{M} corresponding to variables, 
#' 					respectively, to their coefficients
#' @param thresh	(numeric) value used to check whether the influence of the a coefficient
#' 					to reducing the error sum of squares is small enough to conclude that the
#' 					corresponding column in \eqn{X'X} is a linear combination of preceding 
#' 					columns
#' @param tol		(numeric) value used to check numerical equivalence to zero
#' @param Ncpu		(integer) number of cores to be used for parallel processing
#'                  (not yet used)
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @return (list) with two elements:\cr
#' 			\item{SSQ}{(numeric) vector of ANOVA sum of squares}
#' 			\item{LC}{(integer) vector indicating linear dependence of each column}
#' 
#' @references 
#' Goodnight, J.H. (1979), A Tutorial on the SWEEP Operator, The American Statistician, 33:3, 149-158

Csweep <- function(M, asgn, thresh=1e-12, tol=1e-12, Ncpu=1)
{
	nr <- nrow(M)
	stopifnot(nr == ncol(M))
	LC <- ZeroK <- rep(0, length(asgn))
	
	ind <- is.nan(M) | is.na(M) | M == Inf | M == -Inf;
	if(any(ind))
		M[which(ind)] <- 0
	
	swept <- .C("Tsweep", M=as.double(t(M)), k=as.integer(asgn), thresh=as.double(thresh), 
			NumK=as.integer(length(asgn)), nr=as.integer(nr), LC=as.integer(LC), 
			tol=as.double(tol), SSQ=as.double(rep(0, length(unique(asgn)))), 
			PACKAGE="VCA")
	
	res <- list(SSQ=swept$SSQ, 
			LC=tapply(swept$LC, asgn, function(x) length(which(x==1))))
	
	return(res)
}


#' Calling C-implementation of the SWEEP-Operator for Matrix-Inversion
#' 
#' Function calls a fast C-implementation of the SWEEP operator using the
#' transpose of the matrix to be swept.
#' 
#' Transposing prior to applying the SWEEP-operator speeds up things since the 
#' complete matrix is stored in memory in consecutive manner. 
#' This version of the SWEEP-operator is intended for matrix inversion only, thus,
#' not computing ANOVA sum of squares and number of linear dependencies (see function
#' \code{\link{Csweep}}).
#' 
#' This is an utility-function not intended to be called directly.
#' 
#' @param M			(matrix) matrix, representing the augmented matrix \eqn{X'X}
#' @param tol		(numeric) value used to check numerical equivalence to zero
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @return (Matrix) object corresponding to the inverted matrix

Sinv <- function(M, tol=.Machine$double.eps)
{
	nr <- nrow(M)
	stopifnot(nr == ncol(M))
	
#	ind <- is.nan(M) | is.na(M) | M == Inf | M == -Inf;
#	if(any(ind))
#		M[which(ind)] <- 0
#	
	swept <- .C("TsweepFull", M=as.double(t(M)), nr=as.integer(nr), 
			tol=as.double(tol*max(abs(diag(M)))), PACKAGE="VCA")
	
	if(class(M) == "matrix")
		return(matrix(round(swept$M, abs(log10(tol))), nrow=nr, ncol=nr, byrow=TRUE))
	else
		return(Matrix(round(swept$M, abs(log10(tol))), nrow=nr, ncol=nr, byrow=TRUE))
}



#' ANOVA Sum of Squares via Sweeping.
#' 
#' Compute ANOVA Type-1 sum of squares for linear models.
#' 
#' This function performs estimation of ANOVA Type-1 sum of squares
#' using the SWEEP-operator (see reference), operating on the augmented
#' matrix \eqn{X'X}, where \eqn{X} represents the design matrix not differentiating
#' between fixed and random factors.
#' 
#' This is an utility function not intended to be called directly.
#' For each term in the formula the design-matrix \eqn{Z} is constructed.
#' Matrix \eqn{X} corresponds to binding all these \eqn{Z}-matrices together column-wise.
#' 
#' Degrees of freedom for each term are determined by subtracting the number of
#' linearly dependent columns from the total number of column in X asigned to a
#' specific term.
#' 
#' @param Data			(data.frame) with the data
#' @param tobj			(terms) object derived from original formula object
#' @param random		(character) vector, optionally containing information about each
#' 						model term, whether it is random or fixed (only used in mixed models)
#' 
#' @return (list) representing the  with variables:\cr
#' 			\item{aov.tab}{basic ANOVA-table with degrees of freedom (DF), SS and MS}
#' 			\item{Lmat}{(list) with components 'Z' and 'A'}
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @references 
#' 
#' Goodnight, J.H., (1979), A Tutorial on the SWEEP Operator, The American Statistician, 33:3, p.149-158
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' res <- VCA:::getSSQsweep(dataEP05A2_1, terms(y~day/run))
#' str(res)
#' }
getSSQsweep <- function(Data, tobj, random=NULL)
{	
	form     <- formula(tobj)
	resp     <- as.character(form)[2]
	fac      <- attr(tobj, "term.labels")
	int      <- attr(tobj, "intercept") == 1 
	
	N        <- nrow(Data)															
	SS       <- numeric()
	DF       <- numeric()
	Lmat     <- list()                                                      	# compute A-matrices, used in constructing covariance-matrix of VCs
	Lmat$Z   <- list()
	Lmat$Zre <- Matrix(nrow=N, ncol=0)
	y        <- Matrix(Data[,resp], ncol=1)
	
	ord <- attr(tobj, "order")
	vars <- rownames(attr(tobj, "factors"))
	fac <- c(fac, "error")
	DF  <- rep(0, length(fac))
	done <- rep(FALSE, length(fac))					# keep track which DFs have been computed
	names(DF) <- fac
	vars <- vars[-1]
	Nvc  <- length(fac) 							# error not included
	N <- nrow(Data)
	asgn <- NULL
	Lmat <- list()
	Lmat$Z <- list()
	
	if(int)
	{
		Lmat$Zre <- Matrix(1, nrow=N, ncol=1)
		colnames(Lmat$Zre) <- "int"
		asgn     <- 0
	}	
	else
		Lmat$Zre <- Matrix(0, nrow=N, ncol=0)
	
	NumK <- c(rep(0, Nvc-1), N)
	
	for(i in 1:Nvc)                 				# construct design-matrices Z for each term, and matrix X (here Zre)                       
	{
		if(i < Nvc)
		{
			tmpMM <- model.matrix(as.formula(paste(resp, "~", fac[i], "-1", sep="")), Data)
			Lmat$Z[[i]] <- Matrix(tmpMM)
			all0 <- apply(Lmat$Z[[i]], 2, function(x) all(x==0))
			if(any(all0))
			{
				ZeroCol <- which(all0)
				Lmat$Z[[i]] <- Lmat$Z[[i]][,-ZeroCol]
			}	
			if(!is.null(random))
				attr(Lmat$Z[[i]], "type") <- ifelse(fac[i] %in% random, "random", "fixed")
			
			attr(Lmat$Z[[i]], "term") <- fac[i]
			
			NumK[i] <- ncol(Lmat$Z[[i]])
			Lmat$Zre    <- Matrix(cbind(as.matrix(Lmat$Zre), as.matrix(Lmat$Z[[i]])))
			asgn <- c(asgn, rep(i, ncol(Lmat$Z[[i]])))
		}
		else
		{
			Lmat$Z[[Nvc]] <- Diagonal(N)							# include error design matrix
			attr(Lmat$Z[[Nvc]], "type") <- "random"
			attr(Lmat$Z[[Nvc]], "term") <- "error"
		}
	}
	
	attr(Lmat$Zre, "assign") <- asgn
	X <- Lmat$Zre
	Xt <- t(X)
	y  <- Matrix(Data[,resp], ncol=1)
	yt <- t(y)
	
	M <- rbind(	cbind(as.matrix(Xt%*%X), as.matrix(Xt%*%y)), 
			cbind(as.matrix(yt%*%X), as.matrix(yt%*%y)))	
	
	uind <- unique(asgn)							# all factors
	SS <- LC <- NULL
	nr <- nrow(M)
	tmpM <- M
	
	Int <- ifelse(int, 1, 0)
	
	V <- rep(1, ncol(M))
	swept <- Csweep(M, asgn=asgn)					# sweep matrix M
	LC    <- swept$LC
	SS	  <- swept$SSQ
	
	if(int)											# no DF-adjustment for intercept
		LC <- LC[-1]
	
	if(Nvc > 1)
	{
		DF[1:(Nvc-1)] <- NumK[1:(Nvc-1)]-LC				# adjust for linearly dependent variables --> resulting in degrees of freedom
		DF["error"] <- tail(NumK,1)-sum(DF[1:(length(DF)-1)]) - Int
	}
	else
		DF["error"] <- NumK - Int
	
	if(!int)
		SS <- c(M[nr, nr], SS)
	
	SSQ <- abs(diff(SS))
	SSQ <- c(SSQ, tail(SS,1))
	
	aov.tab <- data.frame(DF=DF, SS=SSQ, MS=SSQ/DF)	# basic ANOVA-table
	
	return(list(aov.tab=aov.tab, Lmat=Lmat))
}


#' ANOVA Sum of Squares via Quadratic Forms
#' 
#' Compute ANOVA Type-1 sum of squares for linear models.
#' 
#' This function performs estimation of ANOVA Type-1 sum of squares
#' using an approach of expressing them as quadratic forms in \code{y},
#' the column vector of observations. This is an utility function not
#' intended to be called directly.
#' For each term in the formula the design-matrix \code{Z} and the corresponding
#' \code{A}-matrix 
#' Degrees of freedom for each term are determined calling function \code{\link{anovaDF}}.
#' 
#' @param Data			(data.frame) with the data
#' @param tobj			(terms) object derived from original formula object
#' @param random		(character) vector, optionally containing information about each
#' 						model term, whether it is random or fixed (only used in mixed models)
#' 
#' @return (list) representing the  with variables:\cr
#' 			\item{aov.tab}{basic ANOVA-table with degrees of freedom (DF), SS and MS}
#' 			\item{Lmat}{(list) with components 'Z' and 'A'}
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' res <- VCA:::getSSQqf(dataEP05A2_1, terms(y~day/run))
#' str(res)
#' }

getSSQqf<- function(Data, tobj, random=NULL)
{
	form     <- formula(tobj)
	resp     <- as.character(form)[2]
	fac      <- attr(tobj, "term.labels")
	int      <- attr(tobj, "intercept") == 1 
	
	N        <- nrow(Data)															
	SS       <- numeric()
	DF       <- numeric()
	Lmat     <- list()                                                      	# compute A-matrices, used in constructing covariance-matrix of VCs
	Lmat$Z   <- list()
	asgn     <- NULL
	
	if(int)
	{
		Lmat$Zre <- Matrix(1, nrow=N, ncol=1)
		colnames(Lmat$Zre) <- "int"
		asgn     <- 0
	}
	else
		Lmat$Zre <- Matrix(0, nrow=N, ncol=0)
	
	y        <- Matrix(Data[,resp], ncol=1)
	Lmat$A   <- list()
	
	Nvc <- length(fac) + 1
	
	X1 <- NULL																# X1-matrix for building A-matrices continuously growing 
	
	for(i in 1:Nvc)                                                     	# over variance components (including error)                        
	{
		if(i < Nvc)
		{
			Lmat$Z[[i]] <- Matrix(model.matrix(as.formula(paste(resp, "~", fac[i], "-1", sep="")), Data))
			
			all0 <- apply(Lmat$Z[[i]], 2, function(x) all(x==0))
			if(any(all0))
				Lmat$Z[[i]] <- Lmat$Z[[i]][,-which(all0)]
			asgn <- c(asgn, rep(i, ncol(Lmat$Z[[i]])))
			
			if(i == 1)
			{
				Lmat$A[[i]] <- getAmatrix(X1=Matrix(ifelse(int, 1, 0),ncol=1,nrow=N), X2=Lmat$Z[[i]])		# account for intercept
			}
			else
			{
				X1 <- cbind(X1, as.matrix(Lmat$Z[[i-1]]))
				Lmat$A[[i]] <- getAmatrix(X1=Matrix(X1), X2=Lmat$Z[[i]])
			}
			
			if(!is.null(random))
			{
				attr(Lmat$Z[[i]], "type") <- ifelse(fac[i] %in% random, "random", "fixed")
				attr(Lmat$A[[i]], "type") <- ifelse(fac[i] %in% random, "random", "fixed")
			}			
			attr(Lmat$Z[[i]], "term") <- fac[i]
			attr(Lmat$A[[i]], "term") <- fac[i]
			
			Lmat$Zre    <- Matrix(cbind(as.matrix(Lmat$Zre), as.matrix(Lmat$Z[[i]])))
		}
		else
		{
			Lmat$Z[[i]] <- Diagonal(N)                                		# error
			Lmat$A[[i]] <- getAmatrix(X1=getMM(form, Data))
			attr(Lmat$Z[[i]], "type") <- "random"
			attr(Lmat$A[[i]], "type") <- "random"
			attr(Lmat$Z[[i]], "term") <- "error"
			attr(Lmat$A[[i]], "term") <- "error"
		}
		
		SS <- c(SS, as.numeric(t(y) %*% Lmat$A[[i]] %*% y))					# calculate ANOVA sum of squares
	}
	
	attr(Lmat$Zre, "assign") <- asgn
	
	DF <- anovaDF(form=form, Data=Data, Zmat=Lmat$Z, Amat=Lmat$A)			# determine degrees of freedom
	
	aov.tab <- data.frame(DF=DF, SS=SS, MS=SS/DF)
	rownames(aov.tab) <- c(attr(tobj, "term.labels"), "error")
	
	return(list(aov.tab=aov.tab, Lmat=Lmat))
}



#' ANOVA-Type Estimation of Mixed Models.
#' 
#' Estimate/Predict random effects employing ANOVA-type estimation and obtain generalized least squares estimates
#' of fixed effects for any linear mixed model including random models and linear models.
#' 
#' A Linear Mixed Model, noted in standard matrix notation, can be written as {y = Xb + Zg + e}, where
#' \eqn{y} is the column vector of observations, \eqn{X} and \eqn{Z}{Z} are design matrices assigning fixed (\eqn{b}),
#' respectively, random (\eqn{g}) effects to observations, and \eqn{e} is the column vector of residual errors.
#' Whenever there is an intercept in the model, i.e. the substring "-1" is not part of the model formula, the same
#' restriction as in SAS PROC MIXED is introduced setting the last fixed effect equal to zero. Note, that the results
#' of an linear contrasts are not affected by using an intercept or not, except that constrained fixed effects cannot
#' be part of such contrasts (one could use the intercept estimated instead).
#' 
#' Here, no further restrictions on the type of model are made. One can fit mixed models as well as random models, which
#' constitute a sub-set of mixed models (intercept being the only fixed effect). Variables must be either of type "numeric"
#' or "factor". "character" variables are automatically converted to factors and the response variable has to be numeric, of course. 
#' In case that 'class(Data[,i])' is neither one of these three options, an error is issued. 
#' Even simple linear models can be fitted, i.e. models without a random part (without \eqn{Zg}{Zg}) besides the
#' residual errors. In this case, an Analysis of Variance (ANOVA) table is computed in the same way as done by function 'anova.lm'.
#' 
#' One drawback of using ANOVA-type estimation of random effects is, that random effects are independent, i.e they have
#' zero covariance by definition \eqn{cov(g_i,g_j) = 0}. Another one is that estimated variance components may become negative
#' under certain conditions. The latter situation is addressed by setting negative variance estimates equal to zero and adapting
#' ANOVA mean squares (MS) as \eqn{MS = C * VC}, where \eqn{C} is a coefficient matrix and a function of the design matrix \eqn{[X Z]}
#' and \eqn{VC} is the column-vector of adapted variance components. The Satterthwaite approximation of total degrees of freedom 
#' (DF for total variance) will use adapted \eqn{MS}-values. 
#' 
#' Note, that setting negative VCs equal to zero results in a conservative estimate of the total variance, i.e. it will be larger than
#' the estimate including the negative VC(s). Use parameter 'NegVC=TRUE' to explicitly allow negative variance estimates. 
#' 
#' For further details on ANOVA Type-I estimation methods see \code{\link{anovaVCA}}.
#'
#' @param form				(formula) specifying the linear mixed model (fixed and random part of the model),
#' 							all random terms need to be enclosed by brackets. Any variable not being bracketed
#'                      	will be considered as fixed. Interaction terms containing at least one random factor
#'                      	will automatically be random (Piepho et al. 2003).
#' @param Data				(data.frame) storing all variables referenced in 'form', note that variables can only be
#'                          of type "numeric", "factor" or "character". The latter will be automatically converted to "factor".
#' @param by				(factor, character) variable specifying groups for which the analysis should be performed individually,
#' 							i.e. by-processing
#' @param VarVC.method		(character) string specifying whether to use the algorithm given in Searle et al. (1992) which corresponds to \code{VarVC.method="scm"} or in
#' 							Giesbrecht and Burns (1985) which can be specified via "gb". Method "scm" (Searle, Casella, McCulloch)
#'                      	is the exact algorithm but slower, "gb" (Giesbrecht, Burns) is termed "rough approximation"
#' 							by the authors, but sufficiently exact compared to e.g. SAS PROC MIXED (method=type1) which
#' 							uses the inverse of the Fisher-Information matrix as approximation. For balanced designs all
#'                      	methods give identical results, only in unbalanced designs differences occur. 
#' @param SSQ.method		(character) string specifying the method used for computing ANOVA Type-1 sum of squares and respective degrees of freedom.
#' 							In case of "sweep" funtion \code{\link{getSSQsweep}} will be called, otherwise, function \code{\link{getSSQqf}}
#' @param NegVC         	(logical) FALSE = negative variance component estimates (VC) will be set to 0 and they will not 
#' 							contribute to the total variance (as done e.g. in SAS PROC NESTED, conservative estimate of total variance). 
#' 							The original ANOVA estimates can be found in element 'VCoriginal'. 
#'                      	The degrees of freedom of the total variance are based on adapted mean squares (MS) (see details).
#' 							TRUE = negative variance component estimates will not be set to 0 and they will contribute to the total 
#' 							variance (original definition of the total variance).
#' @param quiet				(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' @return (VCA) object
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @references	
#' 
#' Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York	
#' 
#' Goodnight, J.H. (1979), A Tutorial on the SWEEP Operator, The American Statistician, 33:3, 149-158
#' 
#' Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#' Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486
#' 
#' H.P.Piepho, A.Buechse and K.Emrich (2003), A Hitchhiker's Guide to Mixed Models for Randomized Experiments,
#' J.Agronomy & Crop Science 189, p. 310-322
#' 
#' Gaylor,D.W., Lucas,H.L., Anderson,R.L. (1970), Calculation of Expected Mean Squares by the Abbreviated Doolittle and Square Root Methods., 
#' Biometrics 26 (4): 641-655 
#' 
#' SAS Help and Documentation PROC MIXED, SAS Institute Inc., Cary, NC, USA
#' 
#' @seealso \code{\link{anovaVCA}}
#' 
#' @examples 
#' \dontrun{
#' 
#' data(dataEP05A2_2)
#' 
#' # assuming 'day' as fixed, 'run' as random
#' anovaMM(y~day/(run), dataEP05A2_2)
#' 
#' # assuming both as random leads to same results as
#' # calling anovaVCA
#' anovaMM(y~(day)/(run), dataEP05A2_2)
#' anovaVCA(y~day/run, dataEP05A2_2)
#' 
#' # use different approaches to estimating the covariance of 
#' # variance components (covariance parameters)
#' dat.ub <- dataEP05A2_2[-c(11,12,23,32,40,41,42),]			# get unbalanced data
#' m1.ub <- anovaMM(y~day/(run), dat.ub, SSQ.method="qf", VarVC.method="scm")
#' m2.ub <- anovaMM(y~day/(run), dat.ub, SSQ.method="qf", VarVC.method="gb")		# is faster
#' V1.ub <- round(vcovVC(m1.ub), 12)
#' V2.ub <- round(vcovVC(m2.ub), 12)
#' all(V1.ub == V2.ub)
#' 
#' # make it explicit that "gb" is faster than "scm"
#' # compute variance-covariance matrix of VCs 10-times
#' 
#' system.time(for(i in 1:500) vcovVC(m1.ub))	# "scm"
#' system.time(for(i in 1:500) vcovVC(m2.ub))	# "gb"
#' 
#' 
#' # fit a larger random model
#' data(VCAdata1)
#' fitMM1 <- anovaMM(y~((lot)+(device))/(day)/(run), VCAdata1[VCAdata1$sample==1,])
#' fitMM1
#' # now use function tailored for random models
#' fitRM1 <- anovaVCA(y~(lot+device)/day/run, VCAdata1[VCAdata1$sample==1,])
#' fitRM1
#' 
#' # there are only 3 lots, take 'lot' as fixed 
#' fitMM2 <- anovaMM(y~(lot+(device))/(day)/(run), VCAdata1[VCAdata1$sample==2,])
#' 
#' # the following model definition is equivalent to the one above,
#' # since a single random term in an interaction makes the interaction
#' # random (see the 3rd reference for details on this topic)
#' fitMM3 <- anovaMM(y~(lot+(device))/day/run, VCAdata1[VCAdata1$sample==2,])
#' 
#' # fit same model for each sample using by-processing
#' lst <- anovaMM(y~(lot+(device))/day/run, VCAdata1, by="sample")
#' lst
#' 
#' # fit mixed model originally from 'nlme' package
#'  
#' library(nlme)
#' data(Orthodont)
#' fit.lme <- lme(distance~Sex*I(age-11), random=~I(age-11)|Subject, Orthodont) 
#' 
#' # re-organize data for using 'anovaMM'
#' Ortho <- Orthodont
#' Ortho$age2 <- Ortho$age - 11
#' Ortho$Subject <- factor(as.character(Ortho$Subject))
#' fit.anovaMM1 <- anovaMM(distance~Sex*age2+(Subject)*age2, Ortho)
#' 
#' # use simplified formula avoiding unnecessary terms
#' fit.anovaMM2 <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
#' 
#' # and exclude intercept
#' fit.anovaMM3 <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho)
#' 
#' # compare results
#' fit.lme
#' fit.anovaMM1
#' fit.anovaMM2
#' fit.anovaMM3
#' 
#' # are there a sex-specific differences?
#' cmat <- getL(fit.anovaMM3, c("SexMale-SexFemale", "SexMale:age2-SexFemale:age2")) 
#' cmat
#' 			 
#' test.fixef(fit.anovaMM3, L=cmat)
#' 
#' # former versions of the package used R-function 'lm' and 'anova',
#' # which is significantly slower for sufficiently large/complex models
#' data(realData)
#' datP1 <- realData[realData$PID==1,]
#' system.time(anova.lm.Tab <- anova(lm(y~lot/calibration/day/run, datP1)))
#' # using the quadratic forms approach for estimating ANOVA Type-1 sums of squares
#' system.time(anovaMM.Tab1  <- anovaMM(y~lot/calibration/day/run, datP1, SSQ.method="qf"))
#' # using the sweeping approach for estimating ANOVA Type-1 sums of squares
#' # this is now the default setting (Note: only "gb" method works as VarVC.method)
#' # Also see ?anovaVCA for a comparison of the computational efficiency of "qf" and "sweep".
#' system.time(anovaMM.Tab2  <- anovaMM(y~lot/calibration/day/run, datP1, SSQ.method="sweep"))
#' 
#' # compare results, note that the latter corresponds to a linear model,
#' # i.e. without random effects. Various matrices have already been computed,
#' # e.g. "R", "V" (which are identical in this case).
#' anova.lm.Tab
#' anovaMM.Tab1
#' anovaMM.Tab2
#' }
#' 
#' @seealso \code{\link{anovaVCA}}, \code{\link{VCAinference}}, \code{\link{remlVCA}}, \code{\link{remlMM}}
#' 			\code{\link{ranef}}, \code{\link{fixef}}, \code{\link{vcov}}, \code{\link{vcovVC}}, 
#' 			\code{\link{test.fixef}}, \code{\link{test.lsmeans}}, \code{\link{plotRandVar}}

anovaMM <- function(form, Data, by=NULL, VarVC.method=c("gb", "scm"), SSQ.method=c("sweep", "qf"), 
		NegVC=FALSE, quiet=FALSE)
{
	if(!is.null(by))
	{
		stopifnot(is.character(by))
		stopifnot(by %in% colnames(Data))
		stopifnot(is.factor(by) || is.character(by))
		
		levels  <- unique(Data[,by])
		res <- lapply(levels, function(x) anovaMM(form=form, Data[Data[,by] == x,], NegVC=NegVC, SSQ.method=SSQ.method, VarVC.method=VarVC.method, quiet=quiet))
		names(res) <- paste(by, levels, sep=".")
		return(res)
	}
	
	stopifnot(class(form) == "formula")
	stopifnot(is.data.frame(Data))
	stopifnot(nrow(Data) > 2)                                               	# at least 2 observations for estimating a variance
	
	VarVC.method <- match.arg(VarVC.method)
	SSQ.method   <- match.arg(SSQ.method)
	VarVC.method <- ifelse(SSQ.method == "sweep", "gb", VarVC.method)			# always use "gb", since A-matrices will not be computed	
	
	res <- list()
	res$call <- match.call()
	res$data <- Data															# as provided 
	org.form <- form
	tobj <- terms(form, simplify=TRUE, keep.order=TRUE)                    		# expand nested factors if necessary retain order of terms in the formula
	
	if(length(attr(tobj, "term.labels")) == 0)									# handle pure-error models with 'anovaVCA'
		return(anovaVCA(form, Data))
	
	int   <- res$intercept <- attr(tobj, "intercept") == 1						# has intercept
	form  <- formula(tobj)
	res$terms <- tobj
	
	if(!attr(tobj, "response"))
		stop("You need to include a response variable in the fixed effects formula!")
	resp <- as.character(form)[2]
	res$response <- resp
	
	stopifnot(resp %in% colnames(Data))
	stopifnot(is.numeric(Data[,resp]))
	
	rf <- gregexpr("\\([[:alnum:]]*\\)", as.character(org.form)[3])				# check for random effects
	
	if(rf[[1]][1] != -1)														# identify random variables
	{
		len <- attr(rf[[1]], "match.length")
		pos <- rf[[1]]
		tmp <- NULL
		for(i in 1:length(len))
		{
			tmp <- c(tmp, substr(as.character(org.form)[3], pos[i]+1, pos[i]+len[i]-2))				# remember random factors, exclude brackets
		}
		rf <- tmp
	}
	
	Ndata <- nrow(Data)	
	rmInd <- integer()	
	resp.NA <- is.na(Data[,resp])
	
	if(any(resp.NA))
	{    
		rmInd <- c(rmInd, which(resp.NA))
		if(!quiet)
			warning("There are ", length(which(resp.NA))," missing values for the response variable (obs: ", paste(which(resp.NA), collapse=", "), ")!")
		res$resp.NA <- rmInd
	}    
	
	fac  	<- attr(tobj, "term.labels")

	if(length(fac) > 1)
	{
		rf.ind  <- which(apply(sapply(rf, function(x) regexpr(x, fac)), 1, function(x) any(x>0)))		
	}
	else
	{
		if(length(rf) > 0)
		{
			if(rf == fac)
				rf.ind <- 1
			else
				rf.ind <- numeric(0)
		}
	}
	vars    <- rownames(attr(tobj, "factors"))[-1]	                        # remove response
	Nvc     <- length(fac) + 1 
	
	if(length(rf.ind) > 0)													# at least one random term in 'form'
	{
		res$random <- fac[rf.ind]
		res$fixed  <- fac[-rf.ind]
	}
	else
	{
		res$random <- character(0)											# only fixed effects
		res$fixed  <- fac
	}
	
	res$VCnames <- c(res$random, "error")
	res$Nvc  	<- length(res$VCnames)											# error is one additional VC
	res$Type 	<- if(length(res$fixed) == 0)
				"Random Model"
			else
				"Mixed Model"	
	
	for(i in rev(vars))														# check Data for consistency
	{
		if( any(is.na(Data[,i])))
		{
			NAind <- which(is.na(Data[,i]))
			rmInd <- c(rmInd, NAind)
			if(!quiet)
				warning("Variable '", i,"' has ",length(NAind)," missing values (obs: ", paste(NAind, collapse=", "), ")!" )
		}
		if(!class(Data[,i]) %in% c("numeric", "character", "factor"))
			stop("Variable ",i," is neither one of \"numeric\", \"character\" or \"factor\"!")
		
		if(class(Data[,i]) != "numeric")
		{
			if(!quiet && is.character(Data[,i]))
				warning("Convert variable ", i," from \"charater\" to \"factor\"!")
			
			tmp <- Data[,i]
			tmp.num <- suppressWarnings(as.numeric(as.character(tmp)))		# warnings are very likely here
			
			if(!any(is.na(tmp.num)))										# all levels are numbers
			{
				Width <- max(nchar(as.character(unique(Data[,i])))) + 1
				tmp   <- suppressWarnings(formatC(Data[,i], width=Width))
				uLvl  <- unique(Data[,i])
				Data  <- Data[order(tmp),]
				Data[,i] <- factor(as.character(Data[,i]), levels=uLvl[order(unique(tmp))])
			}		
			else
				Data[,i] <- factor(Data[,i])								# automatically orders
		}
	}
	
	rmInd <- unique(rmInd)
	
	if(length(rmInd) > 0)
		Data <- Data[-rmInd,]												
	
	Data <- na.omit(Data[,rownames(attr(tobj, "factors"))])					# get rid of incomplete observations
	Mean <- mean(Data[,resp], na.rm=TRUE)                                   # mean computed after removing incomplete observations
	Nobs <- N <- nrow(Data)
	y    <- matrix(Data[,resp], ncol=1)										# vector of observations
	
	if(SSQ.method == "qf")
		tmp.res <- getSSQqf(Data, tobj, res$random)
	else
		tmp.res <- getSSQsweep(Data, tobj, res$random)
	
	Lmat    <- tmp.res$Lmat
	aov.tab <- tmp.res$aov.tab
	DF		<- aov.tab[,"DF"]
	SS		<- aov.tab[,"SS"]
	
	rownames(aov.tab) <- c(attr(tobj, "term.labels"), "error")
	
	res$Mean <- Mean
	res$Nobs <- Nobs
	res$aov.org <- aov.tab
	
	rf.ind  <- c(rf.ind, nrow(aov.tab))
	
	C <- getCmatrix(form, Data, aov.tab[,"DF"], "SS", MM=Lmat$Zre)          # compute coefficient matrix C in ss = C * s
	# at this point Zre comprises fixed and random effects
	Ci  <- solve(C)
	C2  <- apply(C, 2, function(x) x/DF)                                    # coefficient matrix for mean squares (MS)
	Ci2 <- solve(C2)
	VC  <- VCorg <- as.matrix( Ci %*% SS)                                   # solve for VCs (p.173)
	VCnam <- rownames(aov.tab)
	VCnam[length(VCnam)] <- "error"
	rownames(VC) <- rownames(VCorg) <- VCnam
	colnames(VC) <- colnames(VCorg) <- "Estimate"
	
	VC  <- VC[rf.ind] 														# only use VC for random terms, inter-fixed effects variation not in scope
	aov.tab <- aov.tab[rf.ind,,drop=F]
	
	rownames(aov.tab)[nrow(aov.tab)] <- "error"
	aov.tab$VC <- VC
	
	res$NegVCmsg <- ""
	res$VCoriginal <- aov.tab[, "VC"]
	
	if(!NegVC)
	{
		IndNegVC <- which(aov.tab[,"VC"] < 0)  
		if(length(IndNegVC) > 0)                                            # there are negative VC
		{
			aov.tab[IndNegVC, "VC"] <- 0                                    # set negative VC to 0
			res$NegVCmsg <- "* VC set to 0"
		}
	} 
	
	totVC  <- sum(aov.tab$VC)
	
	aov.tab <- rbind(total=c(NA, NA, NA, totVC), aov.tab) 	
	aov.tab["total", "DF"] <- SattDF(c(C2[rf.ind, rf.ind] %*% aov.tab[-1, "VC"]), 	# will automatically adapt ANOVA-MS if any VCs were set to 0 
			Ci=Ci2[rf.ind, rf.ind, drop=F], DF=DF[rf.ind])  
	
	suppressWarnings(aov.tab <- cbind(aov.tab, SD=sqrt(aov.tab[,"VC"])))    		# warnings suppressed because sqrt of negative numbers doese not exists
	aov.tab <- cbind(aov.tab, "CV[%]"=aov.tab[,"SD"]*100/Mean)
	aov.tab <- cbind(aov.tab, "%Total"=aov.tab[,"VC"]*100/totVC)
	aov.tab <- aov.tab[,c("DF", "SS", "MS", "VC", "%Total", "SD", "CV[%]")]    
	aov.tab <- as.matrix(aov.tab)
	aov.tab <- apply(aov.tab, 1:2, function(x) ifelse(is.nan(x), NA, x))
	
	res$EstMethod <- "ANOVA"
	
	if(Nobs != Ndata)
		res$Nrm <- Ndata - Nobs                            				# save number of observations that were removed due to missing data
	
	Lmat$y <- y															# column vector of observations
	if(int)
	{
		Lmat$X    <- matrix(1, ncol=1, nrow=nrow(Data))					# design matrix of fixed effects: include intercept --> needs a restriction
		colnames(Lmat$X) <- "int"
		fe.assign <- 0													# '0' indicates intercept
	}
	else
		fe.assign <- NULL
	
	fe <- fac[-rf.ind]

	INT <- ifelse(int, "1", "-1")
	if(length(fe) > 0)
		INT <- paste(INT, "+", sep="")

	fixed.form <- paste(resp, "~", INT, paste(fe, collapse="+"), sep="")
	fixed.form <- as.formula(fixed.form)

	old.opt <- options(contrasts=c("contr.SAS", "contr.poly"))
	
	suppressWarnings({
		X1 	<- model.matrix(fixed.form, data=Data)						# with SAS contrasts but without the columns where restrictions apply
		X10	<- apply(X1, 2, function(x) all(x==0))
		X2 <- model.matrix(	fixed.form, data=Data,						# full model matrix with re-setted contrasts
				contrasts.arg=lapply(Data[,sapply(Data, is.factor), drop=FALSE],
						contrasts, contrasts=FALSE))
		X20	<- apply(X2, 2, function(x) all(x==0))
		
		X2.asgn 	<- attr(X2, "assign")									# keep info		
		
		if(any(X10))
			X1 <- X1[,-which(X10),drop=FALSE]
		
		if(any(X20))
		{
			X2 <- X2[,-which(X20),drop=FALSE]
			X2.asgn 	<- X2.asgn[-which(X20)]
		}				
	})
	
	options(old.opt)													# reset contrasts option
	
	fixed.terms <- terms(fixed.form)
	fe.terms 	<- attr(fixed.terms, "term.labels")
		
	X2[,!colnames(X2) %in% colnames(X1)] <- 0							# transfer restriction from X1 to X2
	if(int)
	{
		colnames(X2)[1] <- "int"
		fe.terms <- c("int", fe.terms)
	}
	
	fe.assign <- X2.asgn
	attr(fe.assign, "terms") <- fe.terms
	
	Lmat$X <- X2

	res$fe.assign 	<- fe.assign										# mapping columns of X to fixed terms in the model formula
	res$fixed.terms <- fixed.terms
	
	Lmat$rf.ind <- rf.ind												# indices of random effects
	Lmat$VCall  <- VCorg												# VC-estimates as if all factors were random
	Lmat$C.SS   <- C
	Lmat$C.MS   <- C2 
	Lmat$Ci.SS  <- Ci
	Lmat$Ci.MS  <- Ci2
	
	res$SSQ.method   <- SSQ.method
	res$VarVC.method <- VarVC.method
	res$aov.tab  <- aov.tab
	res$Matrices <- Lmat
	res$balanced <- if(isBalanced(form, Data)) 
				"balanced"  
			else 
				"unbalanced"
	
	class(res)     <- "VCA"
	
	if(length(Lmat$rf.ind) == 1)										# contains only the residual error
	{
		res$Type <- "Linear Model"
		tmp      <- res$aov.org
		tmp 	 <- cbind(tmp, "F value" = tmp[,"MS"]/tmp[nrow(tmp),"MS"])
		tmp		 <- cbind(tmp, "Pr(>F)" = pf(tmp[,"F value"]/tmp[nrow(tmp), "F value"], tmp[, "DF"], tmp[nrow(tmp), "DF"], lower.tail=FALSE))
		tmp[nrow(tmp), c("F value", "Pr(>F)")] <- NA	
		res$aov.org <- tmp
	}
	
	res <- solveMME(res)
	gc(verbose=FALSE)													# trigger garbage collection
	return(res)
}

#' ANOVA Type-I Degrees of Freedom.
#' 
#' Depending on the type of model, e.g. fully-nested, crossed-nested, etc. algorithms
#' are applied which are believed to be reasonably fast for the respective type of model.
#' 
#' This function is not meant to be called directly. It is invoked by functions \code{\link{anovaVCA}}
#' and \code{\link{anovaMM}}.
#' 
#' Computing the rank of a corresponding \code{A}-matrix, which generates ANOVA sum of squares as
#' quadratic form in \eqn{y} is a general method applicable to all types of models. This usually
#' envokes a singular value decomposition of \eqn{A} which makes it rather slow. 
#' Here, we try to speed up things by differentiating three classes of models, 1) fully-nested 
#' models where DFs are computed as the number of factor-levels minus the sum of higher order terms
#' minus 1, 2) models with only main factors (in this case \code{\link{anova.lm}} is used), 3) models
#' with main factors and interaction terms. Main factors DFs are computed employing function \code{\link{anova.lm}}.
#' DFs for interaction terms are computed by determining the rank of \eqn{A}-matrices. There are certain 
#' designs for which \code{anova.lm} becomes very fast (see examples section of \code{\link{anovaMM}}).
#' 
#' @param form			(formula) object specifying the model for which ANOVA DFs are requested
#' @param Data			(data.frame) with all variables appearing in 'form'
#' @param Zmat			(list) of Z-matrices, representing the design matrices of all model-terms, interpreted 
#' 						as fixed effects (number of columns represent number of unique levels)
#' @param Amat			(list) of A-matrices, generating ANOVA sums of squares as quadratic forms in \eqn{y},
#'                      see \code{\link{anovaVCA}} for details
#' @param tol			(numeric) constant, representing a numeric tolerance used in computing the rank of
#'                      \eqn{A}-matrices
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples
#' \dontrun{
#' # fully-nested design
#' data(realData)
#' datP1 <- realData[realData$PID==1,]
#' system.time(anova.lm.Tab1 <- anova(lm(y~lot/calibration/day/run, datP1)))
#' system.time(anovaMM.Tab1  <- anovaMM(y~lot/calibration/day/run, datP1))
#' anova.lm.Tab1
#' anovaMM.Tab1
#' 
#' # use SSQ.method="qf" (based on quadratic forms)
#' system.time(anovaMM.Tab1.qf  <- anovaMM(y~lot/calibration/day/run, datP1, SSQ.method="qf"))
#' 
#' # compute degrees of freedom
#' VCA:::anovaDF( y~lot/calibration/day/run, datP1,
#' 				  Zmat=anovaMM.Tab1.qf$Matrices$Z,
#' 				  Amat=anovaMM.Tab1.qf$Matrices$A)
#' 
#' # design with only main-factors
#' system.time(anova.lm.Tab2 <- anova(lm(y~lot+calibration+day+run, datP1)))
#' system.time(anovaMM.Tab2  <- anovaMM(y~lot+calibration+day+run, datP1))
#' anova.lm.Tab2
#' anovaMM.Tab2
#' 
#' # use SSQ.method="qf" (based on quadratic forms)
#' system.time(anovaMM.Tab2.qf  <- anovaMM(y~lot+calibration+day+run, datP1, SSQ.method="qf"))
#' 
#' # compute degrees of freedom
#' VCA:::anovaDF( y~lot+calibration+day+run, datP1,
#' 				  Zmat=anovaMM.Tab2.qf$Matrices$Z,
#' 				  Amat=anovaMM.Tab2.qf$Matrices$A)
#' 
#' # design with main-factors and interactions
#' system.time(anova.lm.Tab3 <- anova(lm(y~(lot+calibration)/day/run, datP1)))
#' system.time(anovaMM.Tab3  <- anovaMM( y~(lot+calibration)/day/run, datP1))
#' anova.lm.Tab3
#' anovaMM.Tab3
#' 
#' # use SSQ.method="qf" (based on quadratic forms)
#' system.time(anovaMM.Tab3.qf  <- anovaMM(y~(lot+calibration)/day/run, datP1, SSQ.method="qf"))
#' 
#' # compute degrees of freedom
#' VCA:::anovaDF( y~(lot+calibration)/day/run, datP1,
#' 				  Zmat=anovaMM.Tab3.qf$Matrices$Z,
#' 				  Amat=anovaMM.Tab3.qf$Matrices$A)
#' }

anovaDF <- function(form, Data, Zmat, Amat, tol=1e-8)
{
	form <- terms(form, simplify=TRUE, keep.order=TRUE)
	
	fac  <- attr(form, "term.labels")
	
	if(length(fac) == 0)										# handle pure intercept models
		return(nrow(Data)-1)
	fmat <- attr(form, "factors")[-1,,drop=FALSE]				# remove row corresponding to the response
	csum <- apply(fmat, 2, function(x) length(which(x>0)))
	fn   <- all(diff(csum) == 1)								# fully-nested model? ( i-th term part of (i+1)-th term --> csum[i]+1 == csum[i+1] )
	
	DF <- numeric(length(fac))
	Split <- list()
	
	Nlvl <- sapply(Zmat, ncol)									# number of unique factor levels
	names(Nlvl) <- fac
	
	if(fn)														# very fast DF-algorithm applicable
	{
		for(i in 1:length(fac))
		{
			if(i == 1)
				DF[i] <- Nlvl[i] - ifelse(attr(form, "intercept"), 1, 0)
			else
				DF[i] <- Nlvl[i] - sum(DF[1:(i-1)]) - ifelse(attr(form, "intercept"), 1, 0)
		}
	}
	else														# not fully-nested
	{
		ia <- csum > 1											# interactions? 
		mf <- !ia												# main factors
		
		if(!any(ia))											# only main-factors in the model
		{
			DF <- anova(lm(form, Data))
			DF <- DF[-nrow(DF),"Df"]
		}
		else													# main factors and interactions
		{
			if(any(mf))											# get DFs for all main factors
			{
				resp <- as.character(form)[2]					# response variable
				
				tmp <- anova(lm(paste(resp, "~", paste(fac[which(mf)], collapse="+"), 
										ifelse(attr(form, "intercept"),"", "-1"), sep=""), Data))
				
				DF[which(mf)]  <- tmp[-nrow(tmp), "Df"]
			}
			
			if(any(ia))											# apply computationally expensive DF-algorithm on interactions
			{
				for(i in which(ia))
					DF[i] <- rankMatrix(Amat[[i]], tol=tol)
			}			
		}
	}
	
	DF <- c(DF, nrow(Data) - sum(DF) - ifelse(attr(form, "intercept"), 1, 0))									# residual degrees of freedom
	names(DF) <- c(fac, "error")
	attr(DF, "Nlevel") <- Nlvl
	return(DF)
}


#' Solve Mixed Model Equations.
#' 
#' Function solves the Mixed Model Equations (MMEs) to estimate fixed and random effects.
#' It is for internal use only, thus, not exported.
#' 
#' @param obj			... (VCA) object
#' 
#' @return 	(VCA) object, which has additional elements "RandomEffects" corresponding to the column vector 
#' 		   	of estimated random effects, "FixedEffects" being the column vector of estimated fixed effects. 
#' 			Element "Matrices" has additional elements referring to the elements of the MMEs and element
#' 			"VarFixed" corresponds to the variance-covariance matrix of fixed effects.
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples
#' \dontrun{
#' data(dataEP05A2_1)
#' fit <- anovaVCA(y~day/run, dataEP05A2_1, NegVC=TRUE)
#' fit <- solveMME(fit)
#' ranef(fit)
#' }

solveMME <- function(obj)
{
	stopifnot(class(obj) == "VCA")
	mats   	<- obj$Matrices
	V      	<- mats[["V"]]
	
	if(is.null(V))						# will be the case for ANOVA-type fitted models
	{
		obj  <- getV(obj)
		mats <- obj$Matrices
		V	 <- mats[["V"]]
	}
	Z      	<- mats$Zre
	R 	   	<- mats$R
	G      	<- mats$G
	X      	<- Matrix(mats$X)
	y      	<- mats$y
	
	if(is.null(mats$Vi))
		Vi  <- Solve(V)
	else
		Vi	<- mats$Vi
	
	K 		<- Sinv(t(X) %*% Vi %*% X)		# variance-covariance matrix of fixed effects
	T	   	<- K %*% t(X) %*% Vi
	fixed  	<- T %*% y
	
	mats$Vi <- Vi
	mats$T  <- T
	rownames(fixed) <- colnames(X)
	colnames(fixed) <- "Estimate"
	re		<- obj$RandomEffects
	
	if(is.null(Z))
		re <- NULL
	else
	{
		if(is.null(re))
		{
			re  <- G %*% t(Z) %*% Vi %*% (y - X %*% fixed) 
			rownames(re) <- colnames(Z)
			colnames(re) <- "Estimate"
		}
	}
	obj$RandomEffects <- re
	obj$FixedEffects  <- fixed
	obj$Matrices	  <- mats
	obj$VarFixed      <- K
	return(obj)
}


#' Generic Method for Extracting Random Effects from a Fitted Model.
#' @param object		(object)
#' @param ...			additional parameters
#' @seealso \code{\link{ranef.VCA}}

ranef <- function(object, ...)
	UseMethod("ranef")

#' Extract Random Effects from 'VCA' Object.
#' 
#' Extract random effects and possibly apply a transformation to them (standardization,
#' studentization).
#' 
#' Extracting the 'RandomEffects' element of an 'VCA' object if this exists and applying
#' standardization (mean 0, sd 1) or studentization. For studentized random effects 
#' the i-th random effects is divided by the i-th main diagonal element of matrix \eqn{O = GZ^{T}QZG}{O = GZ'QZG},
#' where \eqn{G} is the covariance-matrix of random effects, \eqn{Z} is a design matrix assigning 
#' random effects to observations and matrix \eqn{Q = V^{-1}(I - H)}{Q = V"(I - H)} (see \code{\link{residuals.VCA}} for further details). 
#' 
#' @param object		(VCA) object where random effects shall be extracted
#' @param term			(character) string specifying a term (factor) for which random effects 
#'                      should be extracted, one can also specify an integer which is interpreted
#'                      as i-th element of 'obj$res.assign$terms'
#' @param mode			(character) string or abbreviation specifying whether "raw" residuals
#'                      should be returned or a transformed version c("student" or "standard")
#' @param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' @param ...			additional parameters
#' 
#' @method ranef VCA
#' @S3method ranef VCA
#' 
#' @references 
#' 
#' Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York	
#' 
#' Laird, N.M., Ware, J.H., 1982. Random effects models for longitudinal data. Biometrics 38, 963-974.
#' 
#' Schuetzenmeister, A. and Piepho, H.P. (2012). Residual analysis of linear mixed models using a simulation approach.
#' Computational Statistics and Data Analysis, 56, 1405-1416
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' fit <- anovaVCA(y~day/run, dataEP05A2_1)
#' ranef(fit)
#' 
#' # get variable-specific random effects (REs)
#' # both extract the same REs
#' ranef(fit, "day")
#' ranef(fit, 1)
#' 
#' # get standardized REs
#' ranef(fit, "day:run", "standard")
#' 
#' # or studentized REs
#' ranef(fit, 2, "stu")
#' }

ranef.VCA <- function(object, term=NULL, mode=c("raw", "student", "standard"), quiet=FALSE, ...)
{
	Call <- match.call()
	
	obj <- object
	
	if(is.list(obj) && class(obj) != "VCA")
	{
		if(!all(sapply(obj, class) == "VCA"))
			stop("Only lists of 'VCA' object are accepted!")
		
		obj.len <- length(obj)
		
		if(!"msgEnv" %in% ls(.GlobalEnv))
			msgEnv <<- new.env(parent=emptyenv())
		
		assign("VCAinference.obj.is.list", TRUE, envir=msgEnv)			# indicate that a list-type object was passed intially
		
		if(is.null(term))
		{
			res <- mapply(	FUN=ranef.VCA, object=obj,
					mode=mode[1], SIMPLIFY=FALSE)
		}
		else
		{
			res <- mapply(	FUN=ranef.VCA, object=obj, term=term,
					mode=mode[1], SIMPLIFY=FALSE)
		}
		names(res) <- names(obj)
		
		if(obj.len == 1)			# mapply returns a list of length 2 in case that length(obj) was equal to 1
			res <- res[1]
		
		rm("VCAinference.obj.is.list", envir=msgEnv)
		
		return(res)
	}	
	
	stopifnot(class(obj) == "VCA")
	mode <- match.arg(mode)
	
	ObjNam  <- deparse(Call$object)
	ObjNam2 <- sub("\\[.*", "", ObjNam)
	
	if(is.null(obj$RandomEffects))
	{
		obj  <- solveMME(obj)
		
		if(length(ObjNam2) == 1 && ObjNam2 %in% names(as.list(.GlobalEnv)))
		{
			expr <- paste(ObjNam, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else
		{
			if( !"VCAinference.obj.is.list" %in% names(as.list(msgEnv)) && !quiet)
				warning("Mixed model equations were solved but results could not be assigned to 'VCA' object!")
		}
	}
	
	if(mode == "student" && is.null(obj$Matrices$Q))
	{
		mats <- obj$Matrices
		X  <- mats$X
		T  <- mats$T
		Vi <- mats$Vi
		mats$H <- H  <- X %*% T
		mats$Q <- Q  <- Vi %*% (diag(nrow(H))-H)
		obj$Matrices <- mats
		
		if(length(ObjNam2) == 1 && ObjNam2 %in% names(as.list(.GlobalEnv)))
		{
			expr <- paste(ObjNam, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else
		{
			if( !"VCAinference.obj.is.list" %in% names(as.list(msgEnv)) && !quiet)
				warning("Matrices 'H' and 'Q' were comuted but could not be assigned to 'VCA' object!")
		}
	}
	
	re <- obj$RandomEffects
	nam <- rownames(re)
	if(!is.null(term) && !is.character(term))
	{
		term <- obj$re.assign$terms[as.integer(term)]
	}
	
	if(mode == "standard")
	{
		tmp <- tapply(re, obj$re.assign$ind, scale)
		re  <- matrix(nrow=nrow(re))
		for(i in 1:length(obj$re.assign$terms))
			re[which(obj$re.assign$ind == i)] <- tmp[[i]]
		
		rownames(re) <- nam
	}
	else if(mode == "student")
	{
		G <- getMat(obj, "G")
		Q <- getMat(obj, "Q")
		Z <- getMat(obj, "Z")
		O <- G %*% t(Z) %*% Q %*% Z %*% G
		re <- re / sqrt(diag(O))
		re <- as.matrix(re)
		rownames(re) <- nam
	}
	
	if(is.null(term) || !term %in% obj$re.assign$terms)
	{
		if(!is.null(term) && !term %in% obj$re.assign$terms && !quiet)
		{
			warning("There is no term in the random part of the formula corresponding to specified 'term'!")
		}
		
		ind <- NULL
		for(i in 1:length(obj$re.assign$terms))
			ind <- c(ind, which(obj$re.assign$ind == i))
		
		re <- re[ind,,drop=FALSE]
		
		attr(re, "mode") <- mode
		attr(re, "term") <- "all"
		
		return(re)
	}
	else
	{
		ind <- which(obj$re.assign$ind == which(obj$re.assign$terms == term)) 
		re  <- re[ind,,drop=F]
		attr(re, "mode") <- mode
		attr(re, "term") <- term
		return(re)
	}
}



#' Generic Method for Extracting Fixed Effects from a Fitted Model.
#' @param object		(object)
#' @param ...			additional parameters
#' @seealso \code{\link{fixef.VCA}}

fixef <- function(object, ...)
	UseMethod("fixef")

#' Extract Fixed Effects from 'VCA' Object.
#' 
#' Conviniently extracting the 'FixedEffects' element of an 'VCA' object. 
#' 
#' The default is to return the fixed effects estimates together with their standard errors.
#' If setting 'type="complex"' or to an abbreviation (e.g. "c") additional inferential statistics
#' on these estimates will be returned, i.e. "t Value", "DF" and respective p-value "Pr > |t|". 
#' One can choose one of three denominator degrees of freedom ('ddfm')-methods.
#' 
#' @param object		(VCA) object where fixed effects shall be extracted
#' @param type			(character) string or partial string, specifying whether
#'                      to return "simple" (reduced) or a rather "complex" (more detailed) 
#'                      information about fixed effects
#' @param ddfm			(character) string specifying the method used for computing the 
#'                      degrees of freedom of the t-statistic. Only used when type="complex".
#' 						Available methods are "contain", "residual", and "satterthwaite".
#' @param tol			(numeric) value representing the numeric tolerance use in comparisons, values
#' 						smaller than 'tol' will be considered equal to 0
#' @param quiet			(logical) TRUE = suppress warning messages, e.g. for non-estimable contrasts
#' @param ...			additional parameters
#' 
#' @method fixef VCA
#' @S3method fixef VCA
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' fit <- anovaVCA(y~day/(run), dataEP05A2_1)
#' fixef(fit)
#' 
#' # for complex models it might take some time computing complex output
#' data(VCAdata1)
#' fit <- anovaMM(y~(lot+device)/(day)/(run), VCAdata1[VCAdata1$sample==2,])
#' fixef(fit, "c")
#' }

fixef.VCA <- function(object, type=c("simple", "complex"), ddfm=c("contain", "residual", "satterthwaite"), 
		tol=1e-12, quiet=FALSE, ...)
{
	Call <- match.call()
	
	obj <- object
	
	if(is.list(obj) && class(obj) != "VCA")
	{
		if(!all(sapply(obj, class) == "VCA"))
			stop("Only lists of 'VCA' object are accepted!")
		
		obj.len <- length(obj)
		
		if(!"msgEnv" %in% ls(.GlobalEnv))
			msgEnv <<- new.env(parent=emptyenv())
		
		assign("VCAinference.obj.is.list", TRUE, envir=msgEnv)			# indicate that a list-type object was passed intially
		
		res <- mapply(	FUN=fixef.VCA, obj=obj, type=type[1],
				ddfm=ddfm[1], tol=tol, quiet=quiet,
				SIMPLIFY=FALSE)
		
		names(res) <- names(obj)
		
		if(obj.len == 1)			# mapply returns a list of length 2 in case that length(obj) was equal to 1
			res <- res[1]
		
		rm("VCAinference.obj.is.list", envir=msgEnv)
		
		return(res)
	}	
	
	stopifnot(class(obj) == "VCA")
	type <- match.arg(type)
	if(length(ddfm) > 1 && type == "complex")
	{
		ddfm <- "satterthwaite"
		if(!quiet)
			warning("Note: 'ddfm' not specified, option \"satterthwaite\" was used!")
	}
	ddfm <- match.arg(ddfm)
	
	vc <- vcov(obj)
	se <- suppressWarnings(sqrt(diag(vc)))
	fe <- obj$FixedEffects
	
	if(is.null(fe))								# solve mixed model equations first
	{
		obj  <- solveMME(obj)
		fe   <- obj$FixedEffects
		nam0 <- deparse(Call$object)
		
		nam1 <- sub("\\[.*", "", nam0)			# remove any index-operators 
		if(length(nam1) == 1 && nam1 %in% names(as.list(.GlobalEnv)))
		{
			expr <- paste(nam0, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else			# warning only if not called on list of VCA-objects
		{
			if( !"VCAinference.obj.is.list" %in% names(as.list(msgEnv)) && !quiet)
				warning("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
		}
	}
	nam <- rownames(fe)
	
	if(type == "simple")
	{		
		fe <- matrix(fe, ncol=1)
		colnames(fe) <- "Estimate"
		fe <- cbind(fe, SE=se)		
		rownames(fe) <- nam
	}
	else
	{
		if(is.null(obj$VarCov))
			obj$VarCov <- vcovVC(obj)
		if(is.null(obj$VarFixed))
			obj$VarFixed <- vcov(obj)
		LC  <- diag(nrow(fe))	
		rownames(LC) <- rownames(fe)
		colnames(LC) <- rownames(fe)
		tst <- test.fixef(obj, LC, ddfm=ddfm, quiet=TRUE)
		tst <- cbind(tst, SE=se)
		fe  <- tst[,c(1,5,2,3,4),drop=F]
		NAs <- is.na(fe[,"Estimate"])
		if(any(NAs))
		{
			fe[which(NAs),"Estimate"] <- 0
			fe[which(NAs), 2:5] <- NA
		}
	}
	
	return(fe)
}




#' Least Squares Means of Fixed Effects.
#' 
#' Computes Least Squares Means (LS Means) of fixed effects for fitted mixed models of class 'VCA'.
#' 
#' Function computes LS Means of fixed effects and their corresponding 
#' standard errors. In case of setting argument 'type' equal to "complex" (or any
#' abbreviation) a \eqn{t}-test is performed on each LS Mean, returning degrees
#' of freedom, t-statistic and corresponding p-values. One can choose from one of three
#' denominator degrees of freedom ('ddfm')-methods.
#' 
#' Actually, function \code{\link{test.fixef}} is called with the "no intercept" 
#' version of the fitted model. The "complex" option is significantly slower for unbalanced
#' designs (see \code{\link{test.fixef}} for details). In case that the 'VarCov' element of
#' the 'VCA' object already exists (calling \code{\link{vcovVC}}), which is the most time 
#' consuming part, results can be obtained in less amount of time.
#' 
#' Standard Errors of LS Means are computed as \eqn{TPT^{T}}{T * P * T'}, where \eqn{T} is the
#' LS Means generating contrast matrix and \eqn{P} is the variance-covariance matrix of
#' fixed effects.
#' 
#' Argument \code{at} can be used to modify the values of covariables when computing LS Means and/or
#' to apply different weighting schemes for (fixed) factor varialbes in the model, e.g. when the prevelance
#' of factor-levels differs from a uniform distribution. Usually, if the weighting scheme is not modified,
#' each factor-level will contribute \eqn{1/N} to the LS Mean, where \eqn{N} corresponds to the number of factor-levels. 
#' 
#' Covariables have to be specified as 'name=value', where value can be a vector of length > 1. 
#' Each value will be evaluated for each row of the original LS Means contrast matrix. 
#' If multiple covariables are specified, the i-th element of covariable 1 will be matched with
#' the i-th element of covariable(s) 2...M, where \eqn{M} is the number of covariables in the model.
#' 
#' To apply a different weighting scheme for factor-variables one has to specify 'factor-name=c(level-name_1=value_1,
#' level-name_2=value_2, ..., level-name_N=value_N)'. The sum of all 'value_i' elements must be equal to 1, otherwise,
#' this factor-variable will be skipped issuing a warning. If any levels 'level-name_i' cannot be found for 
#' factor-variable 'factor-name', this variable will also be skipped and a warning will be issued.
#' See the examples section to get an impression of how this works.
#' 
#' @param obj			(VCA) object having at least one fixed effect
#' @param var			(character) string specifying a fixed effects variable for which
#'                      LS Means should be computed, defaults to all fixed effects, i.e. for
#'         				each level of a fixed effects variable ls means will be computed
#' @param type			(character) "simple" = fast version of computing LS means
#' @param ddfm			(character) string specifying the method used for computing the 
#'                      degrees of freedom of the t-statistic. Only used when type="complex".
#' 						Available methods are "contain", "residual", and "satterthwaite".
#' @param at			(list) where each element corresponds either to a (numeric) covariable or
#' 						to a factor-variable for which the weighting scheme should be adjusted.
#' 						See details section for a thorough description of how argument 'at' works
#' 						and also see the examples.
#' @param contr.mat		(logical) TRUE = the LS Means generating contrast-matrix will be added to the
#' 						result as attribute \code{contrasts}
#' @param quiet			(logical) TRUE = suppress warning messages, e.g. for non-estimable contrasts
#' 
#' @return (matrix) with LS Means of fixed effects and respective standard errors,
#'         in case of 'type="complex"'
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples #
#' \dontrun{
#' data(dataEP05A2_2)
#' fit1 <- anovaMM(y~day/(run), dataEP05A2_2)
#' lsmeans(fit1)
#' lsmeans(fit1,, "complex")
#' 
#' # a more complex model
#' data(VCAdata1)
#' fit2 <- anovaMM(y~(lot+device)/(day)/(run), VCAdata1[VCAdata1$sample==2,])
#' lsmeans(fit2, "lot")
#' lsmeans(fit2, "device", "complex")
#' 
#' # pre-computed 'VarCov' element saves time
#' system.time(lsm1 <- lsmeans(fit2, "device", "complex"))
#' fit2$VarCov <- vcovVC(fit2)
#' system.time(lsm2 <- lsmeans(fit2, "device", "complex"))
#' lsm1
#' lsm2 
#' 
#' # simulate some random data 
#' set.seed(212)
#' id <- rep(1:10,10)
#' x <- rnorm(200)
#' time <- sample(1:5,200,replace=T)
#' y <- rnorm(200)+time
#' snp <- sample(0:1,200,replace=T)
#' dat <- data.frame(id=id,x=x,y=y,time=time,snp=snp)
#' dat$snp <- as.factor(dat$snp)
#' dat$id <- as.factor(dat$id)
#' dat$time <- as.numeric(dat$time)
#' dat$sex <- gl(2, 100, labels=c("Male", "Female"))
#' dat$y <- dat$y + rep(rnorm(2, 5, 1), c(100, 100))
#' 
#' fit3 <- remlMM(y~snp+time+snp:time+sex+(id)+(id):time, dat)
#' 
#' # comute standard LS Means for variable "snp"
#' lsmeans(fit3, var="snp")
#' lsmeans(fit3, var="snp", type="c")    # comprehensive output
#' 
#' # compute LS Means at timepoints 1, 2, 3, 4
#' # Note: original LS Means are always part of the output
#' lsmeans(fit3, var="snp", at=list(time=1:4))
#' 
#' # compute LS Means with different weighting scheme
#' # for factor-variable 'sex'
#' lsmeans(fit3, var="snp", at=list(sex=c(Male=.3, Female=.7)))
#' 
#' # combine covariables at some value and altering the
#' # weighting scheme
#' lsmeans(fit3, var="snp", at=list(time=1:4, sex=c(Male=.3, Female=.7)))
#' 
#' # now with comprehensive output and requesting the
#' # LS Means generating contrast matrix
#' lsmeans(fit3, var="snp", type="complex", contr.mat=TRUE,
#'         at=list(time=1:4, sex=c(Male=.3, Female=.7)))
#' }

lsmeans <- function(obj, var=NULL, type=c("simple", "complex"), ddfm=c("contain", "residual", "satterthwaite"), 
					at=NULL, contr.mat=FALSE, quiet=FALSE)
{
	Call <- match.call()
	
	if(is.list(obj) && class(obj) != "VCA")
	{
		if(!all(sapply(obj, class) == "VCA"))
			stop("Only lists of 'VCA' object are accepted!")
		
		obj.len <- length(obj)
		
		if(is.null(var))
		{
			res <- mapply(	FUN=lsmeans, obj=obj, 
					type=type[1], ddfm=ddfm[1], quiet=quiet,
					SIMPLIFY=FALSE)
		}
		else
		{
			res <- mapply(	FUN=lsmeans, obj=obj, var=var, 
					type=type, ddfm=ddfm, quiet=quiet,
					SIMPLIFY=FALSE)
		}
		
		names(res) <- names(obj)
		
		if(obj.len == 1)			# mapply returns a list of length 2 in case that length(obj) was equal to 1
			res <- res[1]
		
		return(res)
	}	
	
	stopifnot(class(obj) == "VCA")
	stopifnot(obj$Type %in% c("Linear Model", "Mixed Model"))		# won't work for random models
	
	type <- match.arg(type)
	
	if(!is.null(var))												# does this fixed effects variable exist
		stopifnot(var %in% obj$fixed)
	
	if(length(ddfm) > 1 && type == "complex")
	{		
		ddfm <- "satterthwaite"
		if(!quiet)
			warning("Note: 'ddfm' not specified, option \"satterthwaite\" was used!")
	}
	ddfm <- match.arg(ddfm)
	
	if(obj$EstMethod == "REML" && ddfm == "contain" && type=="complex")
	{
		ddfm <- "satterthwaite"
		if(!quiet)
			warning("Note: 'ddfm' set to \"satterthwaite\", currently option \"contain\" does not work for REML-estimation!")
	}
	
	if(obj$Type != "Mixed Model" && !obj$intercept)
	{
		if(!quiet)
			warning("There are no fixed effects for which LS Means could be calculated!")
		return(NA)
	}
	
	T <- lsmMat(obj, var=var, quiet=quiet)					# LS Means generating contrast matrix

	id.mat <- NULL											# being non-NULL means that something was done via 'at'
	
	if(!is.null(at))										# LS Means at some fixed values of covariates and/or with different weighting scheme for fixed factor-variables
	{
		Means <- attr(T, "means")

		dat.class <- attr(T, "var.class")

		nam <- names(at)
		
		T2  <- T
		cnT <- colnames(T)
		fe.terms <- attr(obj$fe.assign, "terms")
		
		num.at <- nam[dat.class[nam] == "numeric"]
		fac.at <- nam[dat.class[nam] == "factor"]

		if(all(is.na(c(num.at, fac.at))) && !quiet)
			warning("Argument 'at' was not correctly specified! Neither covariables nor factor variables could be matched!")
			
		if(length(num.at) == 1 && is.na(num.at))								
			num.at <- character()							# ensure feasibility tests below
		if(length(fac.at) == 1 && is.na(fac.at))
			fac.at <- character()
		
		if(length(num.at) > 0)								# if there are any covariables
		{
			at.mat <- matrix(nrow=max(sapply(at[num.at], length)), ncol=length(num.at))
			colnames(at.mat) <- num.at
			
			for(i in 1:length(num.at))						# fill matrix, could also be user-defined
			{
				if(!num.at[i] %in% cnT)
				{
					if(!quiet)
						warning("There is no fixed term ", paste0("'", nam[i],"'!"),"Possible terms are:", paste(cnT, collapse=", "))
					next
				}
				
				if(length(at[[num.at[i]]]) < nrow(at.mat))
				{
					if(!quiet)
						warning("at[[",num.at[i],"]] does not match the max-length element of 'at', it will be replicated as necessary!")
					
					at.mat[,i] <- rep(at[[num.at[i]]], ceiling(nrow(at.mat)/length(at[[num.at[i]]])))[1:nrow(at.mat)]	# replicate
				}
				else
					at.mat[,i] <- at[[num.at[i]]]
			}
		}
		else
			at.mat <- matrix(nrow=length(fac.at), ncol=0)					# empty matrix

		if(length(fac.at) > 0)										# treat factor-variables, i.e. different weighting scheme
		{
			fac.lst <- vector("list", length(fac.at))
			names(fac.lst) <- fac.at
			
			for(i in 1:length(fac.at))								# over all specified factor variables
			{	
				if(sum(at[[fac.at[i]]]) != 1)
				{
					if(!quiet)
						warning("Sum of all coefficients of factor-variable ", paste0("'", fac.at[i],"'"), " is not equal to 1! It will be skipped!")
					next
				}
				
				skip <- FALSE
				
				tmp.mat <- matrix(nrow=nrow(at.mat), ncol=0)
				
				for(j in 1:length(at[[fac.at[i]]]))					# over all levels of the i-th factor variable
				{
					tmp <- matrix(at[[fac.at[i]]][j], ncol=1, nrow=nrow(at.mat))
					colnames(tmp) <- paste0(fac.at[i], names(at[[fac.at[i]]])[j])
					
					if(!colnames(tmp) %in% colnames(T))
					{
						if(!quiet)
							warning("Factor-level ", paste0("'",colnames(tmp), "'"),
									" of variable ",paste0("'", fac.at[i], "'"),
									" does not correspond to a fixed effect!\n  This element of 'at' will be skipped!")
						skip <- TRUE
					}
					tmp.mat <- cbind(tmp.mat, tmp)
				}
				
				if(skip)
					next
				else
				{
					at.mat <- cbind(at.mat, tmp.mat)
					fac.lst[[i]] <- tmp.mat
				}
			}
		}
		at.mat <- na.omit(at.mat)

		id.mat <- matrix(nrow=nrow(T), ncol=ncol(at.mat))
		cn.id  <- colnames(at.mat)
		colnames(id.mat) <- cn.id

		if(length(num.at) > 0)
		{
			for(i in 1:length(num.at))								# identifies combination of covar-levels
				id.mat[,num.at[i]] <- Means[[num.at[i]]]["mean"]

			tmp.ind <- cn.id[!cn.id %in% num.at]					# those columns not corresponding to covariables
		}
		else
			tmp.ind <- cn.id										# only columns corresponding to factor variable weights
	
		if(length(fac.at) > 0)
			id.mat[,tmp.ind] <- T[,tmp.ind]

		if(ncol(at.mat) > 0)								# only if there is anything to evaluate
		{
			for(i in 1:nrow(at.mat))						# for each combination specified
			{
				T3  <- T
				
				if(length(num.at) > 0)						# if there any covariables
				{
					for(j in 1:length(num.at))				# over covariables
					{
						cnTi <- cnT[grepl(num.at[j], cnT)]					# all relevant columns in matrix T					
						
						for(k in 1:length(cnTi))							# over all columns in which the current covariable appears (main-effect or interaction)
						{
							if(grepl(":", cnTi[k]))							# interaction
							{
								splt <- unlist(strsplit(cnTi[k], ":"))		# split interaction into atomic terms
								
								tmp.val <- 1
								
								for(l in 1:length(splt))							# over all terms in the interaction
								{
									if(splt[l] == num.at[j])						# user-specified level	
										tmp.val <- tmp.val * at.mat[i,num.at[j]]
									else
									{
										if(splt[l] %in% names(Means))						# another numeric covariable -> use mean of this value
											tmp.val <- tmp.val * Means[[splt[l]]]["mean"]
										else												# check with rowname
										{
											if(splt[l] %in% rownames(T3))
											{
												tmp.fac <- rep(0, nrow(T3))
												tmp.fac[grepl(splt[l], rownames(T3))] <- 1		# only that LS Mean affected, which is part of the interaction in column cnTi[k]
												tmp.val <- tmp.val * tmp.fac
											}
											else										# not part of the current interaction
												tmp.val <- tmp.val * rep(0, nrow(T3))
										}
									}
								}
								T3[,cnTi[k]] <- tmp.val
							}
							else
								T3[,cnTi[k]] <- at.mat[i,num.at[j]]					# main effects		
						}					
					}
				}
				
				if(length(fac.at) > 0)
				{
					for(j in 1:length(fac.at))									# over factor variables
					{
						tmp.row <- fac.lst[[fac.at[j]]][i,]
						T3[, tmp.ind] <- rep(tmp.row, rep(nrow(T3), length(tmp.row)))  
					}
				}
				
				T2 <- rbind(T2, T3)
				
				tmp.id.mat <- matrix(rep(at.mat[i,], nrow(T)), ncol=ncol(at.mat), byrow=TRUE)
				colnames(tmp.id.mat) <- colnames(at.mat)
				
				id.mat <- rbind(id.mat, tmp.id.mat)
			}
		}		
		T <- T2										# T2 might be identical to T
	}

	attr(T, "var.class") <- NULL
	attr(T, "means") 	 <- NULL

	x   <- obj
	if(x$intercept)
	{
		x$intercept <- FALSE		
	}
		
	if(type == "complex")						# complex output --> slower
	{	
		if(is.null(obj$VarCov))
			obj$VarCov <- vcovVC(obj)			# determine vcov of VCs if missing
		lsm <- test.fixef(obj, T, ddfm=ddfm, quiet=TRUE, lsmeans=TRUE)
		rownames(lsm) <- rownames(T)
	}
	else										# simple output --> faster
	{
		vc  <- vcov(obj)	
		vc  <- T %*% vc %*% t(T)
		se  <- sqrt(diag(vc))
		lsm <- T %*% obj$FixedEffects
		lsm <- cbind(as.matrix(lsm), SE=se)
	}
	
	if(!is.null(id.mat) && ncol(id.mat) > 0)
		lsm <- cbind(id.mat, lsm)
	
	if(contr.mat)
		attr(lsm, "contrasts") <- T
	
	return(lsm)
}


#' Contrast Matrix for LS Means.
#' 
#' Function determines appropriate contrast matrix for computing the LS Means of
#' each factor level of one or multiple fixed effects variables. This functions implements
#' the 5 rules given in the documentation of SAS PROC GLM for computing the LS Means.
#' 
#' The LS Means correspond to marginal means adjusted for bias introduced by unbalancedness.
#' 
#' @param obj			(VCA) object
#' @param var			(character) string specifyig the fixed effects variable for which
#'                      the LS Means generating matrices should be computed
#' @param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' 
#' @return	(matrix) where each row corresponds to a LS Means generating contrast
#'          for each factor level of one or multiple fixed effects variable(s)
#' 
#' @author Andre Schutzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' fit1 <- anovaMM(y~day/run, dataEP05A2_1)
#' 
#' VCA:::lsmMat(fit1, "day")	# function not exported
#' VCA:::lsmMat(fit1, "run")
#' VCA:::lsmMat(fit1)			# is equal to listing all fixed terms
#' 
#' # a more complex and unbalanced model
#' data(VCAdata1)
#' datS1 <- VCAdata1[VCAdata1$sample == 1, ]
#' set.seed(42)
#' datS1ub <- datS1[-sample(1:nrow(datS1))[1:25],]
#' fit2 <- anovaMM(y~(lot+device)/day/(run), datS1ub)
#' VCA:::lsmMat(fit2, c("lot", "device"))
#' }

lsmMat <- function(obj, var=NULL, quiet=FALSE)
{
	stopifnot(class(obj) == "VCA")
	if(!is.null(var))
		stopifnot(var %in% obj$fixed)
	X <- getMat(obj, "X")
	
	if(is.null(var))
		var <- obj$fixed
	terms 	  <- attr(obj$terms, "factors")
	dat.cls   <- sapply(obj$data[,rownames(terms)[-1]], class)	
	variables <- names(dat.cls)
	fe.assign <- obj$fe.assign						# assignment by integers
	fe.terms  <- attr(fe.assign, "terms")	

	fe.class  <- list()

	for(i in 1:length(fe.terms))
	{
		if(fe.assign[i] == 0)
		{
			fe.class[[i]] <- NULL
		}
		tmp <- unlist(strsplit(fe.terms[i], ":"))
		fe.class[[i]] <- dat.cls[tmp]
	}
	names(fe.class) <- fe.terms

	if(!all(var %in% fe.terms))
		stop("At least one element of 'var' is not a fixed term!")
	
	fe.tvec   <- fe.terms[fe.assign+ifelse(obj$intercept, 1, 0)]				# assignment by names
	
	terms.cls <- rep("factor", length(fe.terms))
	names(terms.cls) <- fe.terms
	fe.tab    <- table(fe.tvec)
	
	lsmm <- matrix(nrow=0, ncol=ncol(X))				
	cn  <- colnames(lsmm) <- colnames(X)
	rn  <- Means <- NULL 
	
	if(any(dat.cls == "numeric"))					# find non-dummy variable columns in X
	{
		num.var <- names(dat.cls[which(dat.cls == "numeric")])
		for(i in 1:length(terms.cls))
		{
			if(any(sapply(num.var, grepl, fe.terms[i])))
			{
				terms.cls[i] <- "numeric"	
			}
		}		
	}
	
	fe.cls 		<- terms.cls[fe.assign+ifelse(obj$intercept, 1, 0)]
	num.terms 	<- terms.cls == "numeric"

	if(any(num.terms))														# are there any numeric terms
	{
		ind <- which(num.terms)
		Means <- list()

		for(i in 1:length(ind))
		{
			tmp		 <- numeric()
			tmp      <- c(tmp, Nlev=as.numeric(fe.tab[fe.terms[ind[i]]]))	# number of columns in X for the current numeric term
			splt     <- unlist(strsplit(fe.terms[ind[i]], ":"))
			tmp.cls  <- dat.cls[splt]
			num.var  <- which(tmp.cls == "numeric")
			fac.var  <- which(tmp.cls == "factor")
			tmp.mean <- apply(obj$data[,names(num.var), drop=FALSE], 1, function(x) eval(parse(text=paste(x, collapse="*"))))	# multiply each value of multiple numeric variables for each row		
			tmp 	 <- c(tmp, mean=mean(tmp.mean))							# mean of current combination of numeric variables
			
			if(length(fac.var) > 0)											# at least one factor variable
			{
				for(j in 1:length(fac.var))									# determine number of levels for each factor variable
				{
					exprs <- paste("c(tmp,",names(fac.var[j]),"=length(unique(obj$data[,\"",names(fac.var[j]),"\"])))", sep="") 
					tmp   <- eval(parse(text=exprs))
				}
			}
			eval(parse(text=paste("Means[[\"", fe.terms[ind[i]], "\"]] <- tmp", sep="")))	# add to list 'Means' and use name of the numeric variable as name of the list element
		}
	}

	for(i in 1:length(var))									# over all fixed terms for which LS Means shall be computed
	{
		tsplt <- unlist(strsplit(var[i], ":"))				# LS Means can only be generated for pure factor variables or combinations of such, no numeric variable may interact
		
		if(any(dat.cls[tsplt] == "numeric"))
		{
			if(!quiet)
				warning("'",var[i],"' is \"numeric\", LS Means cannot be estimated,'",var[i],"' will be skipped!")
			next
		}
		lvl.ind  <- which(fe.tvec == var[i])				# indices of all columns representing levels of var[i]
		lvl.nam  <- cn[lvl.ind]								# use column names of X instead of indices
		lvl.asgn <- unique(fe.assign[lvl.ind])				# value of the fixed effects assignment indicator for the current term var[i]
		tmp.lsm  <- matrix(0, nrow=0, ncol=ncol(X))
		colnames(tmp.lsm) <- cn
		rem.ind.init  <- which(fe.tvec != var[i])			# (init)ial (rem)aining columns which need to be treated differently
		rem.nam.init  <- cn[rem.ind.init]
		rem.asgn.init <- fe.assign[rem.ind.init] 
	
		for(j in 1:length(lvl.ind))							# over levels (columns) of the current factor
		{
			rem.ind  <- rem.ind.init						# re-set for each level of the current (i-th) term
			rem.nam  <- rem.nam.init
			rem.asgn <- rem.asgn.init
			
			con.mat <- matrix(0, nrow=1, ncol=ncol(X))
			colnames(con.mat) <- cn
			splt <- unlist(strsplit(lvl.nam[j], ":"))		# get all sub-terms
			
			if(obj$intercept)
			{
				con.mat[1,"int"] <- 1
				rem.ind  <- rem.ind[ -which(cn == "int")]
				rem.nam  <- rem.nam[ -which(cn == "int")]
				rem.asgn <- rem.asgn[-which(cn == "int")]
			}
			covar <- fe.cls[rem.ind] == "numeric"			# covariates involved?

			if(any(covar))									# [rule 1] (SAS PROC GLM documentation contrast matrix for LS Means)
			{
				cov.ind   <- rem.ind[ which(covar)]			# column-index of current covariate in contrast matrix
				rem.ind   <- rem.ind[-which(covar)]
				rem.nam   <- rem.nam[-which(covar)]
				rem.asgn  <- rem.asgn[-which(covar)]
				cov.asgn  <- fe.assign[cov.ind]				# covariate-indices
				ucov.asgn <- unique(cov.asgn)				

				for(k in 1:length(ucov.asgn))				# over all terms representing covariates
				{
					tmp.term <- fe.terms[ucov.asgn[k] + ifelse(obj$intercept, 1, 0)]
					tmp.info <- Means[[tmp.term]]
					tmp.ind  <- cov.ind[which(cov.asgn == ucov.asgn[k])]
					cn.splt  <- t(sapply(cn[tmp.ind], function(x) unlist(strsplit(x, ":"))))
	
					if(nrow(cn.splt) == 1)													# atomic covariate use mean value
					{
						con.mat[,which(fe.assign == ucov.asgn[k])] <- tmp.info["mean"]
						next
					}
					
					cov.cont <- apply(cn.splt, 1, function(x) any(splt %in% x))				# any sub-terms of the current term found in the levels of the current covariate
					
					if(any(cov.cont))														# covariate is distributed according to a term that is also a sub-term  of the current factor-level 
					{
						con.mat[1,tmp.ind[which(cov.cont)]] <- tmp.info["mean"]/length(which(cov.cont))
					}
					else																	# covariate independent of the current factor-level
					{
						con.mat[1,tmp.ind] <- tmp.info["mean"]/tmp.info["Nlev"]
					}
				}
			}
			
			tmp.splt  <- unlist(strsplit(lvl.nam[j], ":"))									# [rule 2] handling all effects that are contained by the current effect
			contained <- sapply(rem.nam, function(x) 
						all(unlist(strsplit(x, ":")) %in% tmp.splt))				
			
			if(any(contained))
			{
				tmp.asgn  <- rem.asgn[which(contained)]										# this might be multiple elements and this might be different values, e.g. 
				contained <- rem.nam[which(contained)]
				rem.ind   <- rem.ind[-which(rem.asgn %in% tmp.asgn)]						# column belongig to the same term must not be hanled again --> remove them from remaining elements
				rem.nam   <- rem.nam[-which(rem.asgn %in% tmp.asgn)]
				rem.asgn  <- rem.asgn[-which(rem.asgn %in% tmp.asgn)] 
				con.mat[,contained] <- 1
			}
			
			con.mat[1,lvl.ind[j]] <- 1														# [rule 3] setting the columns corresponding to the current effect to 1, all others remain 0
			
			contain <- sapply(rem.nam, function(x) 											# [rule 4] consider effects that contain the current effect
						all(tmp.splt %in% unlist(strsplit(x, ":"))))  
			
			if(any(contain))
			{
				tmp.asgn  <- rem.asgn[which(contain)]
				utmp.asgn <- unique(tmp.asgn)
				
				for(k in 1:length(utmp.asgn))	
				{
					tmp.ind <- rem.ind[which(contain)]
					tmp.ind <- tmp.ind[which(tmp.asgn == utmp.asgn[k])]						# only those indices that belong to the k-th term (assignment value)
					con.mat[1, tmp.ind] <- 1 / length(tmp.ind)								
				}
				
				rem.ind  <- rem.ind[ -which(rem.asgn %in% utmp.asgn)]
				rem.nam  <- rem.nam[ -which(rem.asgn %in% utmp.asgn)]
				rem.asgn <- rem.asgn[-which(rem.asgn %in% utmp.asgn)]
			}
			
			if(length(rem.ind) > 0)															# [rule 5] there are still remaining, not yet treated effects
			{																				#          set these to 1/number of levels
				ufac.asgn <- unique(rem.asgn)
				
				for(k in 1:length(ufac.asgn))
				{
					con.mat[1, which(fe.assign == ufac.asgn[k])] <- 1/fe.tab[fe.terms[ufac.asgn[k] + ifelse(obj$intercept, 1, 0)] ]
				}
			}
			
			tmp.lsm <- rbind(tmp.lsm, con.mat)		
		}
		rownames(tmp.lsm) <- lvl.nam
		lsmm <- rbind(lsmm, tmp.lsm)
	}
	
	attr(lsmm, "var.class") <- dat.cls
	
	if(!is.null(Means))
		attr(lsmm, "means")	<- Means

	return(lsmm)
}




#' Extract Fixed Effects from 'VCA' Object.
#' 
#' For coneniently using objects of class 'VCA' with other packages expecting this
#' function, e.g. the 'multcomp' package for general linear hypotheses for parametric
#' models.
#' 
#' @param object		(VCA) object where fixed effects shall be extracted
#' @param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' @param ...			additional parameters
#'
#' @method coef VCA 
#' @S3method coef VCA
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' fit1 <- anovaMM(y~day/(run), dataEP05A2_1)
#' coef(fit1)
#' fit2 <- anovaVCA(y~day/run, dataEP05A2_1)
#' coef(fit2)
#' }

coef.VCA <- function(object, quiet=FALSE, ...)
{
	Call <- match.call()
	obj <- object
	stopifnot(class(obj) == "VCA")
	fe  <- fixef(obj)
	if(is.null(fe))								# solve mixed model equations first
	{
		obj  <- solveMME(obj)
		fe   <- fixef(obj)
		nam  <- as.character(as.list(Call)$object)
		if(length(nam) == 1 && nam %in% names(as.list(.GlobalEnv)))
		{
			expr <- paste(nam, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else
		{
			if(!quiet)
				warning("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
		}
	}
	fe  <- fe[,"Estimate", drop=F]
	nam <- rownames(fe)
	fe  <- c(fe)
	names(fe) <- nam
	return(fe)
}


#' Calculate Variance-Covariance Matrix of Fixed Effects for an 'VCA' Object.
#' 
#' Return the variance-covariance matrix of fixed effects for a linear mixed model
#' applicable for objects of class 'VCA'.
#' 
#' Actually this function only extracts this matrix or, if not available, calls function
#' \code{\link{vcovFixed}} which performs calculations. It exists for compatibility reasons,
#' i.e. for coneniently using objects of class 'VCA' with other packages expecting this
#' function, e.g. the 'multcomp' package for general linear hypotheses for parametric
#' models.  
#' 
#' @param object		 	(VCA) object for which the variance-covariance matrix of
#'                          fixed effects shall be calculated
#' @param quiet				(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' @param ...				additional parameters
#' 
#' @return (matrix) corresponding to the variance-covariance matrix of fixed effects
#' 
#' @method vcov VCA
#' @S3method vcov VCA
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' fit1 <- anovaMM(y~day/(run), dataEP05A2_1)
#' vcov(fit1)
#' 
#' fit2 <- anovaVCAy~day/run, dataEP05A2_1)
#' vcov(fit2)
#' }

vcov.VCA <- function(object, quiet=FALSE, ...)
{
	obj <- object
	stopifnot(class(obj) == "VCA")
	if(!obj$intercept && length(obj$fixed) == 0)
	{
		if(!quiet)
			warning("There is no variance-convariance matrix of fixed effects for this object!")
		return(NA)
	}
	if(is.null(obj$VarFixed))
		return(vcovFixed(obj))
	else
		return(obj$VarFixed)
}


#' Extract Degrees of Freedom from Linear Hypotheses of Fixed Effects or LS Means.
#' 
#' Determine degrees of freedom for custom linear hypotheses of fixed effects or LS Means
#' using one of three possible approximation methods. 
#' 
#' This is a convenience function to determine the DFs for linear hypotheses in the same way
#' as function \code{\link{test.fixef}}. Only the "DF" part is returned here which can be passed
#' to other functions expecting DFs as input.
#' 
#' @param obj			(VCA) object
#' @param L				(matrix) specifying one or multiple linear hypothese, as returned by function
#'                      \code{\link{getL}} 
#' @param method		(character) the method to be used to determine the degrees of freedom for a 
#'                      linear hypothesis
#' @param ...			additional parameters
#' 
#' @return (numeric) vector with the DFs for each row of 'L'
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' \dontrun{
#' data(VCAdata1)
#' tmpDat <- VCAdata1[VCAdata1$sample==1,]
#' tmpDat <- tmpDat[-c(11,51,73:76),]
#' fitMM <- anovaMM(y~(lot+device)/(day)/(run), tmpDat)
#' fitMM
#' L <- getL(fitMM, c("lot1-lot2", "device1-device2"))
#' getDF(fitMM, L)						# method="contain" is Default
#' getDF(fitMM, L, method="res")
#' 
#' getDF(fitMM, L, method="satt")		# takes quite long for this model
#' }

getDF <- function(obj, L, method=c("contain", "residual", "satterthwaite"), ...)
{
	stopifnot(class(obj) == "VCA")
	method <- match.arg(method)
	return(test.fixef(obj, L=L, ddfm=method, onlyDF=TRUE))
}


#' Calculate Variance-Covariance Matrix and Standard Errors of Fixed Effects for an 'VCA' Object.
#' 
#' The variance-covariance matrix of fixed effects for the linear mixed model in 'obj' is calculated.
#' 
#' The variance-covariance matrix of fixed effects for a linear mixed model corresponds to matrix
#' \eqn{(X^{T}V^{-1}X)^{-}}{(X'V"X)`}, where >\eqn{^{T}}{'}< denotes the transpose operator, >\eqn{^{-1}}{"}< 
#' the regular matrix inverse, and >\eqn{^{-}}{`}< the generalized (Moore-Penrose) inverse of a matrix.
#' 
#' @param obj			(VCA) object for which the variance-covariance matrix of
#'                      fixed effects shall be calculated
#' @param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' 
#' @return (matrix) corresponding to the variance-covariance matrix of fixed effects
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' fit1 <- anovaMM(y~day/(run), dataEP05A2_1)
#' vcov(fit1)
#' 
#' fit2 <- anovaVCA(y~day/run, dataEP05A2_1)
#' vcov(fit2)
#' }

vcovFixed <- function(obj, quiet=FALSE)
{
	Call <- match.call()
	if(!obj$intercept && length(obj$fixed) == 0)
	{
		if(!quiet)
			warning("There is no variance-convariance matrix of fixed effects for this object!")
		return(NA)
	}
	
	X  <- getMat(obj, "X")
	if(is.null(obj$Matrices$Vi))			# MME not yet been solved
	{
		obj  <- solveMME(obj)
		nam  <- as.character(as.list(Call)$obj)
		if(length(nam) == 1 && nam %in% names(as.list(.GlobalEnv)))
		{
			expr <- paste(nam, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else
		{
			if(!quiet)
				warning("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
		}
	}
	Vi <- getMat(obj, "Vi")
	VCov <- obj$VarFixed
	if(is.null(VCov))
		VCov <- Sinv(t(X) %*% Vi %*% X)
	#VCov <- MPinv(t(X) %*% Vi %*% X)
	rownames(VCov) <- colnames(VCov) <- rownames(obj$FixedEffects)
	return(VCov)
}


#' Variance-Covariance Matrix of Fixed Effects as Function of Covariance Parameter Estimates.
#' 
#' This is a helper function for function \code{\link{test.fixef}} approximating degrees of freedom for 
#' linear contrasts of fixed effect parameter estimates.
#' 
#' @param obj			(VCA) object
#' @param x				(numeric) vector of covariance parameter estimates
#' 
#' @return (matrix) corresponding to the variance-covariance matrix of fixed effects

DfSattHelper <- function(obj, x)
{
	stopifnot(class(obj) == "VCA")
	
	rf.ind <- obj$Matrices$rf.ind
#	Zi     <- obj$Matrices$Z
	Zi <- obj$Matrices$Zre
	Z <- Vs <- nam <- NULL 
	
#	for(i in 1:length(rf.ind))									# constructing complete Z-matrix from VC-wise Z-matrices
	
	ura <- unique(obj$re.assign[[1]])
	
	for(i in 1:length(x))
	{
		if(i != length(x))
		{
			ind <- which(obj$re.assign[[1]] == i)
#		Z <- cbind(Z, as.matrix(Zi[[rf.ind[i]]]))
			Z <- cbind(Z, as.matrix(Zi[,ind]))
#		Vs <- c(Vs, rep(x[i], ncol(Zi[[rf.ind[i]]])))
			Vs <- c(Vs, rep(x[i], length(ind)))
#		nam <- c(nam, colnames(Zi[[i]]))
			nam <- c(nam, colnames(Zi)[ind])
		}
		else
		{
			Z  <- cbind(Z, diag(obj$Nobs))
			Vs <- c(Vs, rep(x[i], obj$Nobs)) 
		}
	}
	
	G  <- diag(Vs)												# estimated variance-covariance matrix of random effects
	Z  <- Matrix(Z)
	R  <- getMat(obj, "R")
	X  <- getMat(obj, "X")
	V  <- Z %*% G %*% t(Z) + R
	Vi <- Solve(V)
	P  <- Sinv(t(X) %*% Vi %*% X)
#	P  <- MPinv(t(X) %*% Vi %*% X)
	P <- as.matrix(P)	
	return(P)
}


#' Perform t-Tests for Linear Contrasts on LS Means.
#' 
#' Perform custom hypothesis tests on Least Squares Means (LS Means) of fixed effect.
#' 
#' This function is similar to function \code{\link{test.fixef}} and represents a convenient way of specifying
#' linear contrasts of LS Means. 
#' 
#' @param obj			(VCA) object
#' @param L				(matrix) specifying one or multiple custom hypothesis tests as linear contrasts of LS Means.
#'                      Which LS Means have to be used is inferred from the column names of matrix \eqn{L}, which need to 
#'                      be in line with the naming of LS Means in function \code{\link{lsmeans}}.
#' @param ddfm			(character) string specifying the method used for computing the denominator
#'                      degrees of freedom of t-tests of LS Means. Available methods are "contain", 
#' 						"residual", and "satterthwaite".
#' @param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @seealso \code{\link{test.fixef}}, \code{\link{lsmeans}}
#' 
#' @examples
#' \dontrun{
#' data(dataEP05A2_2)
#' ub.dat <- dataEP05A2_2[-c(11,12,23,32,40,41,42),]
#' fit1 <- anovaMM(y~day/(run), ub.dat)
#' fit2 <- remlMM(y~day/(run), ub.dat)
#' lsm1 <- lsmeans(fit1)
#' lsm2 <- lsmeans(fit2)
#' lsm1
#' lsm2
#' 
#' lc.mat <- getL(fit1, c("day1-day2", "day3-day6"), "lsm")
#' lc.mat[1,c(1,2)] <- c(1,-1)
#' lc.mat[2,c(3,6)] <- c(1,-1)
#' lc.mat
#' test.lsmeans(fit1, lc.mat) 
#' test.lsmeans(fit2, lc.mat)
#' 
#' # fit mixed model from the 'nlme' package
#'  
#' library(nlme)
#' data(Orthodont)
#' fit.lme <- lme(distance~Sex*I(age-11), random=~I(age-11)|Subject, Orthodont) 
#' 
#' # re-organize data for using 'anovaMM'
#' Ortho <- Orthodont
#' Ortho$age2 <- Ortho$age - 11
#' Ortho$Subject <- factor(as.character(Ortho$Subject))
#' 
#' # model without intercept
#' fit.anovaMM <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho)
#' fit.remlMM1 <- remlMM( distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho)
#' fit.remlMM2 <- remlMM( distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho, cov=FALSE)
#' lsm0 <- lsmeans(fit.anovaMM)
#' lsm1 <- lsmeans(fit.remlMM1)
#' lsm2 <- lsmeans(fit.remlMM2)
#' lsm0
#' lsm1
#' lsm2
#' 
#' lc.mat <- matrix(c(1,-1), nrow=1, dimnames=list("int.Male-int.Female", c("SexMale", "SexFemale")))
#' lc.mat
#' test.lsmeans(fit.anovaMM, lc.mat)	
#' test.lsmeans(fit.remlMM1, lc.mat)
#' test.lsmeans(fit.remlMM2, lc.mat)
#' }	

test.lsmeans <- function(obj, L, ddfm=c("contain", "residual", "satterthwaite"), quiet=FALSE)
{
	stopifnot(class(obj) == "VCA")
	lsmm <- lsmMat(obj, NULL, quiet=quiet)
	stopifnot(all(colnames(L) %in% rownames(lsmm)))			# only existing LS Means can be used
	
	lsmm <- lsmm[which(colnames(L) %in% rownames(lsmm)),]
	lsmm <- as.matrix(lsmm)
	
	U <- matrix(nrow=0, ncol=ncol(lsmm))
	colnames(U) <- colnames(lsmm)
	
	for(i in 1:nrow(L))
	{
		U <- rbind(U, matrix(L[i,], nrow=1) %*% lsmm)
	}
	
	if(is.null(obj$VarCov))
		obj$VarCov 	<- vcovVC(obj)
	
	res <- test.fixef(obj, U, ddfm=ddfm, quiet=quiet)
	rownames(res) <- rownames(L)
	return(res)
}


#' Perform t-Tests for Linear Contrasts on Fixed Effects.
#' 
#' This function performs t-Tests for one or multiple linear combinations (contrasts) of estimated 
#' fixed effects.
#' 
#' Here, the same procedure as in SAS PROC MIXED ddfm=satterthwaite (sat) is implemented. 
#' This implementation was inspired by the code of function 'calcSatterth' of R-package 'lmerTest'. 
#' Thanks to the authors for this nice implementation. \cr
#' Note, that approximated Satterthwaite degrees of freedom might differ from 'lmerTest' and SAS PROC MIXED.
#' Both use the inverse Fisher-information matrix as approximation of the variance-covariance matrix
#' of variance components (covariance parameters). Here, either the exact algorithm for ANOVA-estimators of
#' variance components, described in Searle et. al (1992) p. 176, or the approximation presented in Giesbrecht and 
#' Burns (19985) are used. For balanced designs their will be no differences, usually. 
#' In case of balanced designs, the Satterthwaite approximation is equal to the degrees of freedom of the highest
#' order random term in the model (see examples).
#' 
#' @param obj			(VCA) object
#' @param L				(numeric) vector or matrix, specifying linear combinations of the fixed effects, in the latter case,
#'                      each line represents a disctinct linear contrast
#' @param ddfm			(character) string specifying the method used for computing the denominator
#'                      degrees of freedom for tests of fixed effects or LS Means. Available methods are
#' 						"contain", "residual", and "satterthwaite".
#' @param method.grad	(character) string specifying the method to be used for approximating the gradient
#'                      of the variance-covariance matrix of fixed effects at the estimated covariance parameter
#'                      estimates (see function 'grad' (numDeriv) for details)
#' @param tol			(numeric) value specifying the numeric tolerance for testing equality to zero
#' @param quiet			(logical) TRUE = suppress warning messages, e.g. for non-estimable contrasts
#' @param opt			(logical) TRUE = tries to optimize computation time by avoiding unnecessary computations
#'                      for balanced datasets (see details). 
#' @param onlyDF		(logical) TRUE = only the specified type of degrees of freedom are determined without carrying out
#' 						the actual hypothesis test(s)
#' @param ...			further parameters (for internal use actually)
#'  
#' @return (numeric) vector or matrix with 4 elements/columns corresponding to "Estimate", "t Value", "DF", and
#'         "Pr > |t|".
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com} inspired by authors of R-package 'lmerTest'
#' 
#' @references 
#' 
#' Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York	
#' 
#' Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#' Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486
#' 
#' SAS Help and Documentation PROC MIXED (MODEL-statement, Option 'ddfm'), SAS Institute Inc., Cary, NC, USA
#' 
#' @seealso \code{\link{test.lsmeans}}, \code{\link{getL}}
#' 
#' @examples
#' \dontrun{
#' data(dataEP05A2_2)
#' ub.dat <- dataEP05A2_2[-c(11,12,23,32,40,41,42),]
#' fit1 <- anovaMM(y~day/(run), ub.dat)
#' fit2 <- remlMM(y~day/(run), ub.dat)
#' fe1 <- fixef(fit1)
#' fe1
#' fe2 <- fixef(fit2)
#' fe2
#' lc.mat <- getL( fit1, c("day1-day2", "day3-day6"))
#' lc.mat
#' test.fixef(fit1, lc.mat, ddfm="satt") 
#' test.fixef(fit2, lc.mat, ddfm="satt")
#' 
#' # some inferential statistics about fixed effects estimates
#' L <- diag(nrow(fe1))
#' rownames(L) <- colnames(L) <- rownames(fe1)
#' test.fixef(fit1, L)
#' test.fixef(fit2, L)
#' 
#' # using different "residual" method determining DFs
#' test.fixef(fit1, L, ddfm="res")
#' test.fixef(fit2, L, ddfm="res")  
#' 
#' # having 'opt=TRUE' is a good idea to save time 
#' # (in case of balanced designs)
#' data(VCAdata1)
#' datS3 <- VCAdata1[VCAdata1$sample==3,]
#' fit3 <- anovaMM(y~(lot+device)/(day)/run, datS3)
#' fit4 <- remlMM(y~(lot+device)/(day)/run, datS3)  
#' fit3$VarCov <- vcovVC(fit3)
#' fe3 <- fixef(fit3)
#' fe4 <- fixef(fit4)
#' L <- diag(nrow(fe3))
#' rownames(L) <- colnames(L) <- rownames(fe)
#' system.time(tst1 <- test.fixef(fit3, L))
#' system.time(tst2 <- test.fixef(fit3, L, opt=FALSE))
#' system.time(tst3 <- test.fixef(fit4, L, opt=FALSE))
#' tst1
#' tst2
#' tst3
#' }

test.fixef <- function(	obj, L, ddfm=c("contain", "residual", "satterthwaite"),
		method.grad="simple", tol=1e-12, quiet=FALSE, opt=TRUE,
		onlyDF=FALSE, ...)
{
	stopifnot(class(obj) == "VCA")
	ddfm <- match.arg(ddfm)
	
	if(obj$EstMethod == "REML" && ddfm == "contain")
	{
		ddfm <- "satterthwaite"
		if(!quiet)
			warning("Note: 'ddfm' set to \"satterthwaite\", currently option \"contain\" does not work for REML-estimation!")
	}
	
	args <- list(...)
	lsmeans <- args$lsmeans						# this function is called when LS Means with complex output is requested
	if(is.null(lsmeans))
		lsmeans <- FALSE
	
	b <- fixef(obj)[,"Estimate", drop=F]
	
	if(is.null(dim(L)))
		L <- matrix(L, nrow=1)
	
	stopifnot(ncol(L) == nrow(b))
	
	if(nrow(L) > 1)
	{
		res <- matrix(ncol=5, nrow=nrow(L))
		colnames(res) <- c("Estimate", "DF", "SE", "t Value", "Pr > |t|")
		for(i in 1:nrow(L))
			res[i,] <- test.fixef(obj=obj, L=L[i,,drop=F], ddfm=ddfm, method.grad=method.grad, quiet=quiet, opt=opt, onlyDF=onlyDF, lsmeans=lsmeans)
		rownames(res) <- rownames(L)
		if(onlyDF)
			return(res[,"DF"])
		attr(res, "ddfm") <- ddfm
		return(res)
	}	
	else
	{
		ind <- which(L != 0)						# test whether fixef(obj, "complex") was called
		
		if(length(ind) == 1 && any(abs(b[ind,]) < tol) && !lsmeans)
		{
			if(!quiet)
				warning("Linear contrast'", L ,"' not estimable!")
			res <- rep(NA, 4)
			names(res) <- c("Estimate", "DF", "t Value", "Pr > |t|")
			return(res)
		}		
		est		<- as.numeric(L %*% b)
		sgn		<- sign(est)
		X		<- getMat(obj, "X")
		P    	<- vcov(obj)
		lPl  	<- L %*% P %*% t(L)
		lPli 	<- solve(lPl)
		r       <- as.numeric(rankMatrix(lPli))
		SVD 	<- eigen(lPl)		
		items   <- list(SVD=SVD, r=r, b=b)	
				
		se <- L %*% P %*% t(L)						# compute standard error of linear contrast
		se <- sqrt(diag(se))

		DF		<- getDDFM(obj, L, ddfm, tol=tol, method.grad=method.grad, opt=opt, items=items)
		if(onlyDF)
		{
			return(c(NA, DF, NA, NA))
		}
		t.stat 	<- as.numeric(sqrt((t(L %*% b) %*% lPli %*% (L %*% b))/r))
		
		res <- matrix(c(est, DF, se, sgn*t.stat, 2*pt(abs(t.stat), df=DF, lower.tail=FALSE)), nrow=1,
				dimnames=list(NULL, c("Estimate", "DF", "SE", "t Value", "Pr > |t|"))) 
		
		attr(res, "ddfm") <- ddfm
		
		return( res )
	}	
}


#' Construct Linear Contrast Matrix for Hypothesis Tests.
#' 
#' Function constructs coefficient/contrast matrices from a string-representation of linear hypotheses.
#' 
#' Function constructs matrices expressing custom linear hypotheses of fixed effects or
#' LS Means. The user has to specify a string denoting this contrast which is then 
#' transformed into a coefficient/contrast matrix. This string may contain names of fixed effects
#' belonging to same same fixed term, numeric coefficients and mathematical operators "+"
#' and "-" (see examples).
#' 
#' @param obj			(VCA) object
#' @param s				(character) string or vector of strings, denoting one or multiple
#'                      linear contrasts
#' @param what			(character) string specifying whether to construct contrast matrices
#'                      of fixed effects ("fixed") or LS Means ("lsmeans"), abbreviations are allowed.
#' 
#' @return (matrix) representing one linear hypothesis of fixed effects or LS Means per row 
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples
#' \dontrun{
#' data(dataEP05A2_2)
#' fit <- anovaMM(y~day/(run), dataEP05A2_2)
#' L <- getL(fit, c("day1-day2", "day5-day10"), what="fixef")
#' L
#' test.fixef(fit, L=L)
#' 
#' # another custom hypothesis
#' L2 <- getL(fit, "0.25*day1+0.25*day2+0.5*day3-0.5*day4-0.5*day5")
#' 
#' # more complex model
#' data(VCAdata1)
#' dataS2 <- VCAdata1[VCAdata1$sample==2,]
#' fit.S2 <- anovaMM(y~(lot+device)/day/(run), dataS2)
#' L3 <- getL(fit.S2, c("lot1-lot2", "lot1:device3:day19-lot1:device3:day20", 
#' 						"lot1:device1:day1-lot1:device1:day2"))
#' test.fixef(fit.S2, L3)
#' }

getL <- function(obj, s, what=c("fixef", "lsmeans"))
{
	what <- match.arg(what)
	nwhat <- c(fixef="fixed effects", lsmeans="LS Means")
	
	if(what == "fixef")
		b <- fixef(obj, quiet=TRUE)
	else
		b <- lsmeans(obj, quiet=TRUE)
	n <- rownames(b)
	
	if(length(s) > 1)
	{
		L <- try(t(sapply(s, function(x) getL(obj, x, what=what))), silent=TRUE)
		if(class(L) == "try-error")
			stop(L[1])
		colnames(L) <- n
		return(L)
	}	
	contr <- rep(0, nrow(b))
	cname <- names(s)
	s <- gsub("\\*", "", s)
	
	splt <- sapply(unlist(strsplit(s, "\\+")), strsplit, "\\-")
	
	fac <- NULL
	sgn <- NULL
	
	for(i in 1:length(splt))
	{
		if(i == 1)
		{
			if(splt[[1]][1] == "")
				sgn <- -1
			else
				sgn <- 1 
		}
		else
			sgn <- 1									# 1st element always "+" since it was firstly used as split-char
		
		sgn <- c(sgn, rep(-1, length(splt[[i]])-1))
		
		for(j in 1:length(splt[[i]]))
		{
			tmp <- regexpr("^[[:digit:]]*\\.?[[:digit:]]*", splt[[i]][j])
			
			if(tmp > -1 && attr(tmp, "match.length") > 0)				# there is a factor at the beginning of the string
			{
				tmp.fac <- substr(splt[[i]][j], 1, attr(tmp, "match.length"))
				fac <- c(fac, as.numeric(tmp.fac) * sgn[j])	
				splt[[i]][j] <- sub(tmp.fac, "", splt[[i]][j])			# remove factor sub-string
			}
			else
				fac <- c(fac, 1 * sgn[j])
		}
	}
	
	splt <- unlist(splt)
	names(splt) <- NULL
	
	if(!all(splt %in% n))
		stop("\nError: There are terms which do not belong to the '",nwhat[what],"'!")
	
	res <- sapply(n, regexpr, splt)
	rownames(res) <- splt
	cnl <- nchar(colnames(res))	
	
	tms <- apply(res, 1, function(x){
				ind <- which(x > 0)
				len <- cnl[ind]
				return(ind[which(len == max(len))])
			})
	contr[tms] <- fac	
	return(matrix(contr, nrow=1, dimnames=list(cname, n)))
}

#' Degrees of Freedom for Testing Linear Contrast of Fixed Effects and Least Square Means.
#' 
#' There are three methods implemented, which are all available in SAS PROC MIXED, namely 
#' "contain", "residual", and "satterthwaite" approximations. See the documentation of SAS
#' PROC MIXED for details on this topic.
#' 
#' The implementation of the Satterthwaite approximation was inspired by the code of function 
#' 'calcSatterth' of R-package 'lmerTest'.
#' 
#' @param obj			(VCA) object
#' @param L				(numeric) vector specifying the linear combination of the fixed effect or
#'                      LS Means
#' @param ddfm			(character) string specifying the method used for computing the denominator
#'                      degrees of freedom for tests of fixed effects or LS Means. Available methods are
#' 						"contain", "residual", and "satterthwaite".
#' @param tol			(numeric) value specifying the numeric tolerance for testing equality to zero
#' @param method.grad	(character) string specifying the method to be used for approximating the gradient
#'                      of the variance-covariance matrix of fixed effects at the estimated covariance parameter
#'                      estimates (see function 'grad' (numDeriv) for details)
#' @param opt			(logical) TRUE = tries to optimize computation time by avoiding unnecessary computations
#'                      for balanced datasets (see \code{\link{test.fixef}}). 
#' @param items			(list) of pre-computed values
#' 
#' @return (numeric) vector with the specified type of degrees of freedom
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @seealso \code{\link{test.fixef}}

getDDFM <- function(obj, L, ddfm=c("contain", "residual", "satterthwaite"), tol=1e-12, method.grad="simple", opt=TRUE, items=NULL)
{
	stopifnot(class(obj) == "VCA")
	stopifnot(!is.null(colnames(L)))
	ddfm <- match.arg(ddfm)
	
	if(ddfm == "residual")
	{
		return(obj$Nobs - rankMatrix(getMat(obj, "X")))			# as described in documentation of SAS PROC MIXED
	}
	else if(ddfm == "contain")
	{		
		cn <- colnames(L)										# fe or LS Means names
		cn <- cn[which(L[1,] != 0)]
		
		if(length(cn) == 1 && cn == "int")						# handle intercept fixed effect
		{
			return(DF <- min(obj$aov.org[obj$random, "DF"]))
		}
		fe <- obj$fixed
		tmp <- sapply(fe, function(x) gregexpr(x, cn))			# can names of fixed terms be found in column names of L?
		
		if(class(tmp) == "list")								# there was only a single non-zero column in L --> list returned
			fe <- sapply(tmp, function(x) x != -1)
		else													# mulitple non-zero columns in L
		{
			fe <- apply(tmp, 2, function(y) any(y==1))
		}		
		
		if(any(fe))
		{
			fe <- fe[which(fe)]
			fe <- names(fe)
			rn <- obj$random
			rn <- sapply(rn, function(x) all(sapply(fe, grepl, x)))
			if(any(rn))
			{
				DF <- min(obj$aov.org[names(rn), "DF"])
				return(DF)
			}
			else
			{
				tmpZ <-  getMat(obj, "Z")
				
				if(!is.null(tmpZ))
					return(obj$Nobs - rankMatrix(cbind(as.matrix(getMat(obj, "X")), as.matrix(tmpZ))))
				else
					return(obj$Nobs - rankMatrix(as.matrix(getMat(obj, "X"))))
			}
		}
		else
		{
			return(obj$Nobs - rankMatrix(cbind(as.matrix(getMat(obj, "X")), as.matrix(getMat(obj, "Z")))))	# no random term including the fixed effects term --> residual DFs
		}
	}
	else if(ddfm == "satterthwaite")
	{
		stopifnot(class(obj) == "VCA")		
		if(is.null(dim(L)))
			L <- matrix(L, nrow=1)		
		
		if(obj$balanced == "balanced" && opt && obj$EstMethod != "REML" && FALSE)
		{
			return(obj$aov.tab[2,"DF"])			
		}
		r       <- items$r
		SVD 	<- items$SVD
		nu.m 	<- NULL											# see SAS-help for SAS PROC MIXED -> model /ddfm=sat
		A <- obj$VarCov											# same names as in SAS Help of PROC MIXED, Model statement, option "ddfm=sat"
		if(is.null(A))
			A <- vcovVC(obj)
		
		VCs <- obj$aov.tab[-1, "VC"]							# all variance components but total
		
		for( m in 1:length(SVD$values) )
		{   
			g <- grad(function(x)  L %*% DfSattHelper(obj, x) %*% t(L),  VCs, method = method.grad)
			nu.m <- c(nu.m, 2 * (SVD$values[m])^2/(t(g) %*% A %*% g))
		}
		
		E  <- sum( (nu.m/(nu.m-2)) * as.numeric(nu.m>2))
		DF <- 2*E*as.numeric(E>r)/(E-r)
		
		return(DF)
	}
}






#' Calculate Variance-Covariance Matrix of Variance Components of 'VCA' objects.
#' 
#' This function computes the variance-covariance matrix of variance components (VC) either
#' applying the approach given in the \eqn{1^{st}}{1st} reference ('method="scm"') or using
#' the approximation given in the \eqn{2^{nd}}{2nd} reference ('method="gb"').
#' 
#' When 'method="scm"' is used function \code{\link{getVCvar}} is called implementing this rather
#' time-consuming algorithm. Both approaches, respectively the results they generate, diverge for
#' increasing degree of unbalancedness. For balanced designs, they seem to differ only due to 
#' numerical reasons (error propagation).
#' 
#' This function is called on a 'VCA' object, which can be the sole argument. In this case the value
#' assigned to element 'VarVC.method' of the 'VCA' object will be used
#' (see \code{\link{getVCvar}} for computational details).
#' 
#' @param obj			(VCA) object
#' @param method		(character) string, optionally specifying whether to use the algorithm given in the
#' 						1st reference ("scm") or in the 2nd refernce ("gb"). If not not supplied, the 
#' 						option is used coming with the 'VCA' object.
#' @param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' 
#' @return (matrix) corresponding to variance-covariance matrix of variance components
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @references
#' Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York
#' 
#' Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#' Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486
#' 
#' @examples
#' \dontrun{
#' data(realData)
#' dat1 <- realData[realData$PID==1,]
#' fit  <- anovaVCA(y~lot/calibration/day/run, dat1, SSQ.method="qf") 
#' vcovVC(fit)
#' vcovVC(fit, "scm")		# Searle-Casella-McCulloch method (1st reference)
#' vcovVC(fit, "gb")		# Giesbrecht and Burns method (2nd reference)
#' }

vcovVC <- function(obj, method=NULL, quiet=FALSE)
{
	Call <- match.call()
	
	stopifnot(class(obj) == "VCA")
	if(!is.null(method))
		method <- match.arg(method, c("scm", "gb"))
	else
		method <- obj$VarVC.method
	
	VCvar <- obj$VarCov
	if(!is.null(VCvar))
		return(VCvar)
	
	Z  <- obj$Matrices$Z
	
	if(method == "scm")								
	{
		Ci <- getMat(obj, "Ci.SS")
		A  <- obj$Matrices$A
		if(obj$Type %in% c("Linear Model", "Mixed Model"))
			VC <- obj$Matrices$VCall
		else
		{
			VC <- obj$aov.tab[-1, "VC", drop=FALSE]
			VCnam <- rownames(obj$aov.tab)[-1]
			VCnam[length(VCnam)] <- "error"
			rownames(VC) <- VCnam
		}
		VCvar <- getVCvar(Ci=Ci, A=A, Z=Z, VC=VC)
		
		if(obj$Type %in% c("Linear Model", "Mixed Model"))
		{
			ind <- obj$Matrices$rf.ind
			VCvar <- VCvar[ind, ind]
		}
	}
	else
	{
		if(is.null(obj$VarFixed))
		{
			obj  <- solveMME(obj)
			nam  <- as.character(as.list(Call)$obj)
			if(length(nam) == 1 && nam %in% names(as.list(.GlobalEnv)))
			{
				expr <- paste(nam, "<<- obj")		# update object missing MME results
				eval(parse(text=expr))
			}
			else
			{
				if(quiet)
					warning("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
			}
		}		
		VCvar <- getGB(obj)							# apply Giesbrecht & Burns approximation
	}
	attr(VCvar, "method") <- method					# which method was used
	
	return(VCvar)
}


#' Extract a Specific Matrix from a 'VCA' Object.
#' 
#' For convinience only, extracting a specific matrix from the 
#' "Matrices" element of a 'VCA' object if this matrix exists.
#' 
#' When 'mat="Z"' the design matrix of random effects will be returned.
#' If one is interested in the design matrix of random effects for a specific
#' variance component use a name like "Z" + NAME, where NAME has to be equal to
#' the name of the VC in the 'VCA' object (see examples). The same applies to 
#' the A-matrices in the quadratic forms, use "A" + NAME for extracting a specific 
#' A-matrix.
#' 
#' @param obj			... (VCA) object
#' @param mat			... (character) string specifying the matrix to be extracted
#' 
#' @return (matrix) as requested by the user
#' 
#' @examples
#' \dontrun{
#' data(dataEP05A2_1)
#' fit <- anovaVCA(y~day/run, dataEP05A2_1)
#' getMat(fit, "Z")
#' getMat(fit, "Zday")
#' getMat(fit, "Zday:run")
#' getMat(fit, "Zerror")
#' getMat(fit, "V")			 	# Var(y)
#' getMat(fit, "G")				# Var(re)
#' }

getMat <- function(obj, mat)
{
	stopifnot(class(obj) == "VCA")
	mats <- obj$Matrices
	if(is.null(mats))
		return(NULL)
	else
	{		
		if(grepl("^A", mat) || grepl("^Z", mat))
		{
			if(mat == "Z")
			{
				return(mats$Zre)
			}
			else
			{
				for(i in 1:length(mats$Z))					# same length as mats$A
				{
					Znam <- paste("Z", attr(mats$Z[[i]], "term"), sep="")
					Anam <- paste("A", attr(mats$A[[i]], "term"), sep="")
					
					if(mat == Znam)
					{
						return(mats$Z[[i]])
					}
					if(mat == Anam)
					{
						return(mats$A[[i]])
					}
				}
				return(NULL)
			}
		}
		else
		{
			return(mats[[mat]])
		}
	}
}


#' Determine V-Matrix for a 'VCA' Object.
#' 
#' Determine the estimated variance-covariance matrix of observations \eqn{y}.
#' 
#' A linear mixed model can be written as \eqn{y = Xb + Zg + e}, where \eqn{y} is the column
#' vector of observations, \eqn{X} and \eqn{Z} are design matrices assigning fixed (\eqn{b}),
#' respectively, random (\eqn{g}) effects to observations, and \eqn{e} is the column vector of 
#' residual errors.
#' The variance-covariance matrix of \eqn{y} is equal to \eqn{Var(y) = ZGZ^{-T} + R}{Var(y) = ZGZ' + R}, where \eqn{R}
#' is the variance-covariance matrix of \eqn{e} and \eqn{G} is the variance-covariance matrix of \eqn{g}.
#' Here, \eqn{G} is assumed to be a diagonal matrix, i.e. all random effects \eqn{g} are mutually independent
#' (uncorrelated).
#' 
#' @param obj			(VCA) object
#' @return (VCA) object with additional elements in the 'Matrices' element, including matrix \eqn{V}.
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

getV <- function(obj)
{
	stopifnot(class(obj) == "VCA")
	mats <- obj$Matrices
	R <- Diagonal(obj$Nobs) * obj$aov.tab["error", "VC"]
	
	if(as.character(obj$terms)[3] == "1")
	{
		mats$V <- mats$R <- R
		mats$G <- mats$Zre <- NULL
		obj$Matrices <- mats
		return(obj)
	}
	Z <- Vs <- nam <- VCnam <- assign <- NULL
	Zi <- mats$Z
	VC <- mats$VCall									# all VCs, i.e. as if the model were a random model
	rf.ind <- mats$rf.ind
	
	rf.ind <- rf.ind[-length(rf.ind)]					# remove index to error VC
	ind <- 1
	for(i in rf.ind)									# constructing complete Z-matrix from VC-wise Z-matrices
	{
		if(ind == 1)
			Z <- cbind(Z, as.matrix(Zi[[i]]))
		else
			Z <- cbind(as.matrix(Z), as.matrix(Zi[[i]]))
		Vs <- c(Vs, rep(VC[i], ncol(Zi[[i]])))
		nam <- c(nam, colnames(Zi[[i]]))
		VCnam <- c(VCnam , attr(Zi[[i]], "term"))
		assign <- c(assign, rep(ind, ncol(Zi[[i]])))
		ind <- ind + 1									# ind vector specifies which REs belong to which terms in the formula
	}
	
	if(is.null(Z))
	{
		V <- R
		G <- NULL
	}
	else
	{
		Z <- Matrix(Z)
		assign <- list(ind=assign, terms=VCnam)
		G <- Diagonal(ncol(Z)) * Vs						# estimated variance-covariance matrix of random effects
		rownames(G) <- colnames(G) <- nam
		V <- Z %*% G %*% t(Z) + R						# compute V
	}
	colnames(Z) <- nam
	mats$Zre <- Z										# complete Z-matrix
	mats$G	 <- G
	mats$V	 <- V
	mats$R   <- R
	obj$Matrices <- mats
	obj$re.assign <- assign
	
	return(obj)
}



#' Coefficient Matrix for (V)ariance (C)omponent (A)nalysis.
#' 
#' Function \code{getCmatrix} computes the coefficient matrix used in equating observed ANOVA sum of squares (or mean squares) to their expected values
#' as linear combination of the unobservable, true variance components. This "can be viewed in some sense as a special form of the method of moments"
#' approach (Searle et al. 1992, "Variance Components", p. 173).
#' 
#' Functions implements the algorithm for finding coefficient-matrix \eqn{C} of a method of moments approach to ANOVA-estimation
#' of variance components (VC), given in the first reference. Matrix \eqn{C} corresponds to the coefficient matrix equating
#' expected ANOVA mean squares (MS) to observed values as linear combinations of the unknown VCs to be estimated. 
#' The most computationally expensive parts of the algorithm were implemented in \code{C}, speeding up things a significantly.
#' 
#' Consider formulas: \eqn{m_{MS} = Cs}{ms = C * s} and \eqn{m_{SS} = Ds}{ss = D * s}, where \eqn{m_{MS}}{ms} and
#' \eqn{m_{SS}}{ss} are a column-vectors of observed ANOVA MS-, respectively, ANOVA SS-values, \eqn{C} and \eqn{D} are
#' coefficient matrices, equating \eqn{m_{MS}}{ms}, respectively, \eqn{m_{SS}}{ss} to linear combinations of VCs \eqn{s},
#' which are to be estimated. Once matrix \eqn{C} or \eqn{D} is found, pre-multiplying \eqn{s} with its inverse gives ANOVA-estimators of VCs.
#' This function implements the algorithm described in the first reference, the "Abbreviated Doolittle and
#' Square Root Methods". One can convert matrices \eqn{C} and \eqn{D} into the other by multiplying/dividing each element by the respective degrees
#' of freedom (DF) associated with the corresponding factor in the model. If \eqn{\diamond}{~} denotes the operator for element-wise multiplication of two 
#' matrices of the same order (Hadamard-product), \eqn{d} is the column vector of DFs, and \eqn{M} is a \eqn{r \times r}{(r x r)} quadratic matrix with 
#' \eqn{M = d1_r^{T}}{M = d 1_r'}, \eqn{1_r^{T}}{1_r'} being the transpose of a column-vector of ones with \eqn{r} elements, 
#' then holds: \eqn{M = D \diamond M}{C = D ~ M}. 
#' 
#' @param form      (formula) object specifying the random model
#' @param Data      (data.frame) containing all variables in 'form'
#' @param DF        (numeric) vector with the degrees of freedom for each VC, if not
#'                  specified it will be determined automatically by refitting 'form' 
#'                  to 'Data'
#' @param type      (character) "MS" = mean squares coefficient matrix 
#'                              "SS" = sum of squares coefficient matrix
#' @param digits    (integer) numeric tolerance expressed in number of digits. This is used
#'                  for testing whether a value is equal to zero (round(x,digits) == 0).
#' @param MM		(Matrix) object referring to the overparameterized model matrix of the full
#'                  model, if provided, it does not need to be computed twice
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @seealso \code{\link{anovaVCA}}, \code{\link{getMM}}
#' 
#' @references 
#'  Gaylor,D.W., Lucas,H.L., Anderson,R.L. (1970), Calculation of Expected Mean Squares by the Abbreviated Doolittle and Square Root Methods. Biometrics 26(4):641-655
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' C_ms <- getCmatrix(y~day/run, dataEP05A2_1, type="MS")
#' C_ms
#' C_ss <- getCmatrix(y~day/run, dataEP05A2_1, type="SS")
#' C_ss
#' aov.tab <- anova(lm(y~day+day:run, dataEP05A2_1))
#' aov.tab
#' apply(C_ss, 2, function(x) x/aov.tab[,"Df"])
#' }

getCmatrix <- function(form, Data, DF=NULL, type=c("MS", "SS"), digits=8L, MM=NULL)
{
	stopifnot(class(form) == "formula")
	stopifnot(is.data.frame(Data))
	if(as.character(form)[3] == "1")
		return(matrix(DF))
	
	type <- match.arg(type)    
	form <- formula(terms(form, simplify=TRUE, keep.order=TRUE))
	
	if(is.null(DF))
	{
		DF <- anova(lm(form, Data))[,"Df"]
	}
	
	if(is.null(MM))
		MM   <- getMM(form, Data)
	
	asgn <- attr(MM, "assign")
	NVC  <- max(asgn)                           # number of variance components (without mean level)              
	nc   <- ncol(MM)                            # number of effects without intercept
	nr   <- nrow(MM)
	
	amat <- Amat <- Bmat <- matrix(0, nc-1, nc-1)
	
	tol  <- eval(parse(text=paste("1e-", digits, sep="")))
	
	if(inherits(MM, "CsparseMatrix"))			# even faster implementation for sparse matrices
	{	
		res <- .C(	"getAmatBmatSparse", xEl=as.double(MM@x), iEl=as.integer(MM@i),
				pEl=as.integer(MM@p), NobsCol=as.integer(diff(MM@p)), 
				ncol=as.integer(nc), nrow=as.integer(nr), tol=as.double(tol), 
				Amat=as.double(Amat), Bmat=as.double(Bmat), amat=as.double(amat),
				Cmat=double(NVC^2), Nvc=as.integer(NVC), asgn=as.integer(asgn[-1]),
				DF=as.integer(DF), PACKAGE="VCA")
	}
	else										# fast C-implementation for dense matrices
	{
		res <- .C("getAmatBmat", mm=as.double(MM), 										# call C-implementation of the most elaborate part of the algorithm
				ncol=as.integer(nc), nrow=as.integer(nr), tol=as.double(tol), 
				Amat=as.double(Amat), Bmat=as.double(Bmat), amat=as.double(amat),
				Cmat=double(NVC^2), NVC=as.integer(NVC), asgn=as.integer(asgn[-1]),
				DF=as.integer(DF),	PACKAGE="VCA")
	}
	
	C <- matrix(res$Cmat, NVC, NVC)
	C <- rbind(C, rep(0, ncol(C)))					# account for error
	C <- cbind(C, rep(1,nrow(C)))
	
	if(type=="SS")                                  # adapt to sum of squares (SS)
	{
		C <- apply(C, 2, function(x) x <- x*DF)
	}
	
	return(C)
}

#' Determine A-Matrix.
#' 
#' Function \code{getAmatrix} finds A-matrices used to express ANOVA sum of squares as quadratic forms in \eqn{y}{y}.
#' 
#' Compute A-Matrix for expressing ANOVA Type-1 sequential sum of squares \eqn{m_{SS}}{SS} as quadratic form in \eqn{y} (column vector of observations),
#' \eqn{y^{T}Ay = m_{SS}}{y' * A * y = SS}.
#' In random models, when the intercept is the only fixed effect, the design matrix of fixed effects is equal to a column vector
#' where each element is equal to one (1). In this setting, it holds \eqn{1^{T}A1 = 0}{1' * A * 1 = 0}, i.e. all elements of A summing to zero.
#' 
#' This function is not intended to be used directly.
#' 
#' @param X1        (matrix) lower order model matrix (comprising all factors entering model before 'X2'-factor)
#' @param X2        (matrix) higher order model matrix (entering model after 'X1' factors)
#' 
#' @return (matrix) for obtaining quadratic forms in \eqn{y} expression ANOVA sum of squares
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

getAmatrix <- function(X1, X2)
{
	I <- Diagonal(nrow(X1))   
	if(missing(X2))                     					# A-matrix for error
	{
		A <- I-X1 %*%MPinv(t(X1) %*% X1) %*% t(X1) 
	}
	else
	{
		stopifnot(nrow(X1) == nrow(X2))       
		X12 <- Matrix(cbind(as.matrix(X1), as.matrix(X2)))
		
		A1 <- (I-X1  %*% MPinv(t(X1)  %*% X1)  %*% t(X1))
		A2 <- (I-X12 %*% MPinv(t(X12) %*% X12) %*% t(X12)) 
		A <- A1 - A2
	}    
	return(A)
}


#' Moore-Penrose Inverse of a Matrix.
#' 
#' This function is originally defined in package 'MASS'. It was adapted
#' to be able to deal with matrices from the 'Matrix' package, e.g. sparse
#' matrices.
#' 
#' @param X			(object) two-dimensional, for which a Moore-Penrose inverse
#'                  has to be computed
#' @param tol		(numeric) tolerance value to be used in comparisons
#' 
#' @return (object) A Moore-Penrose inverse of X.
#' 
#' @author Authors of the 'MASS' package.

MPinv <- function (X, tol = sqrt(.Machine$double.eps)) 
{
	if (length(dim(X)) > 2L) 
		stop("'X' must be two-dimensional")
	Xsvd <- svd(X)
	Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
	if (all(Positive)) 
		Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
	else if (!any(Positive)) 
		array(0, dim(X)[2L:1L])
	else 
		Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
}



#' Covariance-Matrix of Variance Components.
#' 
#' Function \code{getVCvar} determines the covariance matrix of ANOVA-type estimates of variance components (VC)
#' according to the method given in the reference.
#' 
#' This function implements the (exact) method for computing the variance-covariance matrix of variance components
#' obtained emplyoing ANOVA-type estimation of unbalance data, described in Searle et al. (1992) "Variance Components", Wiley, p.176.
#' One feature of this method is that the asymptotic covariance matrix of VCs produced by SAS PROC MIXED (method=type1)
#' (inverse of the Fisher-Information matrix) is equal to the one computed here, in case of balanced designs (data). 
#' For unbalanced designs, both matrices are likely to differ.
#' 
#' It is for internal use only, thus, not exported.
#' 
#' @param Ci        (matrix) inverted C-matrix of coefficients equating observed Sum of Squares (SS) to expected values.
#' @param A         (list) of A-matrices representing quadratic forms of ANOVA-Type I sums of squares
#' @param Z         (list) of Z-matrices, the design matrices assigning random effects to observations for each variance component.
#' @param VC        (numeric) vector of variance components, i.e. sigma^2.
#' 
#' @return (matrix) covariance matrix of estimated variance components
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @references 
#'                  Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_3)
#' res <- anovaVCA(y~day/run, dataEP05A2_3, SSQ.method="qf")
#' res
#' mat <- res$Matrices
#' Var <- VCA:::getVCvar(Ci=mat$Ci.SS, A=mat$A, Z=mat$Z, VC=res$aov.tab[-1,"VC"])
#' round(Var, 12)
#' }

getVCvar <- function(Ci, A, Z, VC)
{
	Ci <- as.matrix(Ci)
	W  <- matrix(ncol=ncol(Ci), nrow=nrow(Ci))
	VC <- matrix(VC, ncol=1, dimnames=list(rownames(VC), NULL))    
	Zs  <- vector("list", length=length(Z))
	
	for(i in 1:length(Z))
	{
		Zs[[i]] <- Z[[i]] %*% t(Z[[i]])      # compute ZZ' for all Z_i only once
	} 
	
	getWij <- function(Ai, Aj, Zs, VC)       # compute elements of W-matrix
	{        
		AiZ <- vector("list", length=length(Z))
		AjZ <- vector("list", length=length(Z))
		
		for(i in 1:length(Z))
		{
			AiZ[[i]] <- Ai %*% Zs[[i]]
			AjZ[[i]] <- Aj %*% Zs[[i]]
		}
		
		M <- matrix(ncol=length(Z), nrow=length(Z))
		
		for(i in 1:length(Z))
		{
			for(j in i:length(Z))
			{
				M[i,j] <- M[j,i] <- Trace(AiZ[[i]] %*% AjZ[[j]])
			}
		}    
		return(t(VC) %*% M %*% VC)
	}
	
	for(i in 1:nrow(W))                     # compute W-matrix (matrix in outer brackets in 21b on p.176)
	{
		for(j in i:ncol(W))         
		{
			W[i,j] <- W[j,i] <- getWij(Ai=A[[i]], Aj=A[[j]], Zs=Zs, VC=VC)
		}
	}
	
	vc <- 2 * Ci %*% W %*% t(Ci)
	rownames(vc) <- colnames(vc) <- rownames(VC)
	
	return(vc)          # general formula (21b) on p.176, with W being the matrix in outer brackets
}


#' Overparameterized Design Matrices.
#' 
#' Function \code{getMM} constructs overparameterized design matrices from a model formula and a data.frame.
#' 
#' This function constructs the overparameterized design matrix for a given dataset 'Data' according to
#' the model formula 'form'. Each combination of factor-levels and or numeric variables is identified
#' and accounted for by a separate column. See examples for differences compared to function 'model.matrix' (stats).
#' This type of design matrix is used e.g. in constructing A-matrices of quadratic forms in \eqn{y} expressing
#' ANOVA sums of squares as such. This is key functionality of functions \code{\link{anovaVCA}} and \code{\link{anovaMM}}. 
#' 
#' @param form      (formula) with or without response specifying the model to be fit
#' @param Data      (data.frame) with the data
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' \dontrun{
#' # load example data (CLSI EP05-A2 Within-Lab Precision Experiment)
#' data(dataEP05A2_3)
#' tmpData <- dataEP05A2_3[1:10,] 
#' 
#' # check out the differences
#' getMM(~day+day:run, tmpData)
#' model.matrix(~day+day:run, tmpData)
#' 
#' # adapt factor variables in 'tmpData'
#' tmpData$day <- factor(tmpData$day)
#' 
#' # check out the differences now
#' getMM(~day+day:run, tmpData)
#' model.matrix(~day+day:run, tmpData)
#' 
#' # numeric covariate 'cov'
#' tmpData2 <- dataEP05A2_3[1:10,] 
#' tmpData2$cov <- 10+rnorm(10,,3)
#' model.matrix(~day*cov, tmpData2)
#' }
#' 
#' @seealso \code{\link{getAmatrix}}

getMM <- function(form, Data)
{
	stopifnot(class(form) == "formula")
	stopifnot(is.data.frame(Data))
	tform <- terms(form, simplify=TRUE, keep.order=TRUE)
	int   <- attr(tform, "intercept") == 1
	form  <- as.character(tform)
	fmat  <- attr(tform, "factors")
	if(length(fmat) == 0)											# case y~1 (y is response here)
		return(Matrix(1, ncol=1, nrow=nrow(Data)))
	tlabs <- rownames(fmat)[which(apply(fmat, 1, sum) > 0)]			# excludes response if available
	tlab.cls <- sapply(Data[,tlabs,drop=F], class)
	fac.obs  <- tlab.cls[tlab.cls == "factor"]						# factor variables in terms of form
	num.obs  <- tlab.cls[tlab.cls != "factor"]						# non-factor variables
	N <- nrow(Data)
	lvls <- strsplit(form[ifelse(length(form) == 3, 3, 2)], "\\+")  # obtain single factors
	lvls <- gsub(" ", "", unlist(lvls))  
	
	if(int)
	{
		mm <- matrix(1, nrow=N, ncol=1)                             # include intercept
		colnames(mm) <- "int"
		assign <- 0
	}
	else
	{
		mm <- matrix(nrow=N, ncol=0)
		assign <- NULL
	}
	
	for(i in 1:length(lvls))                                        # over terms in the formula
	{   
		if(grepl("-.?1", lvls[i]))
			lvls[i] <- sub("-.?1", "", lvls[i])
		
		if(grepl(":", lvls[i]))                                     # crossed terms
		{
			tmpVar  <- unlist(strsplit(lvls[i], ":"))
			facVar  <- tmpVar[tmpVar %in% names(fac.obs)]
			numVar  <- tmpVar[tmpVar %in% names(num.obs)]
			
			VarNum  <- length(numVar) > 0
			if(length(facVar) > 0)
			{
				tmpData <- Data[,facVar,drop=F]                                         # cols of variables in lvls[i]
				VarFac  <- TRUE
				if(ncol(tmpData) > 1)
				{
					tmpName <- t(apply(tmpData, 1, function(x) paste(facVar, x, sep="")))   # terms = factor levels concatenated to variable names 
					tmpName <- apply(tmpName, 1, paste, collapse=":")                       # terms connected by ":"
					tmpName <- unique(tmpName)                              
				}
				else
				{
					tmpName <- paste(facVar, unique(tmpData[,1]), sep=":")
				}
				fac  <- apply(tmpData, 1, paste, collapse="")	
				Data <- cbind(Data, as.factor(fac))                                     # add new factor-variable to Data
				colnames(Data)[ncol(Data)] <- lvls[i]
			}
		}    
		else														# main factors (no interaction)
		{
			if(lvls[i] %in% names(fac.obs))
			{
				tmpName <- paste(lvls[i], unique(Data[,lvls[i]]), sep="")   
				VarNum  <- FALSE
				VarFac  <- TRUE
			}
			else
			{
				tmpName <- numVar <- lvls[i]
				VarNum  <- TRUE
				VarFac  <- FALSE
			}
		}
		
		if(VarFac)
		{
			eff <- unique(Data[,lvls[i]])
			Ni  <- length(eff)                                          # number of effects for the i-th factor
			tmp <- matrix(0, N, Ni)
			
			colnames(tmp) <- unique(tmpName)
			
			for(j in 1:Ni)                                              # over effects of the i-th factor 
			{
				tmp[which(Data[,lvls[i]] == eff[j]),j] <- 1
			}
		}
		else
		{
			Ni  <- 1
			tmp <- matrix(1,nrow=N)
			colnames(tmp) <- tmpName
		}        
		
		if(VarNum)										
		{
			for(j in 1:ncol(tmp))
			{
				tmp[,j] <- tmp[,j] * Data[,numVar]	
			}
		}
		mm <- cbind(mm, tmp)
		assign <- c(assign, rep(i, Ni))
	}
	mm <- Matrix(mm)
	attr(mm, "assign") <- assign
	
	return(mm)
}



#' Satterthwaite Approximation for Total Degrees of Freedom and for Single Variance Components.
#' 
#' This function estimates degrees of freedom of the total variance (type="total")
#' in random models or individual variance components (type="individual"). 
#' It bases on the results of the unified approach to ANOVA-type estimation 
#' of variance components as implemented in functions \code{\link{anovaVCA}} 
#' and \code{\link{anovaMM}}.
#' 
#' Function is used internally, thus, it is not exported. Option 'type="total"' is used in 
#' functions \code{\link{anovaVCA}} and \code{\link{anovaMM}} for approximating total DF.
#' Option 'type="individual"' is used in function \code{\link{VCAinference}} when choosing
#' 'ci.method="satterthwaite"' for approximating DFs for individual variance components.
#' 
#' @param MS        	(numeric) vector of sequential mean squares (ANOVA type-1).
#' @param Ci       		(matrix) where elements are numeric values representing the inverse of the coefficient
#'                      matrix for calculation of expected mean squares (see \code{\link{anovaVCA}}, \code{\link{getCmatrix}}).
#' @param DF        	(numeric) vector with the degrees of freedom for each factor in a ANOVA type-1 model.
#' @param type			(character) string specifying whether "total" degrees of freedom should be approximated or those of
#' 						individual variance components
#'  
#' @return numeric value representing the Satterthwaite DFs of the total variance.
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' 
#' \dontrun{
#' data(dataEP05A2_2)
#' res <- anovaVCA(y~day/run, dataEP05A2_2)
#' VCA:::SattDF(res$aov.tab[-1,"MS"], getMat(res, "Ci.MS"), res$aov.tab[-1,"DF"], type="tot")
#' 
#' # now approximating individual DF for variance components
#' VCA:::SattDF(res$aov.tab[-1,"MS"], getMat(res, "Ci.MS"), res$aov.tab[-1,"DF"], type="i")
#' }

SattDF <- function(MS, Ci, DF, type=c("total", "individual"))
{
	type <- match.arg(type)
	if(type == "individual")
	{
		res <- NULL
		for(i in 1:nrow(Ci))
			res <- c(res, SattDF(MS, Ci[i,,drop=FALSE], DF))
		return(res)
	}
	I <- matrix(1,nrow(Ci),1)
	t <- t(I) %*% Ci
	t <- c(t) * c(MS)
	
	r1 <- sum(t)^2 / sum(t^2/DF)        # DF total
	
	return(r1)  
}


#' Standard Printing Method for Objects of Class 'VCA'.
#' 
#' Function prints 'VCA' objects as returned by function \code{\link{anovaVCA}}.
#' 
#' @param x         (VCA) object of class 'VCA' as returned by function 'anovaVCA'.
#' @param digits    (integer) number of digits numeric values are rounded to before printing.
#' @param ...       additional arguments to be passed to or from methods.
#' 
#' @method print VCA
#' @S3method print VCA

print.VCA <- function(x, digits=6L, ...)
{
	EstMeth <- x$EstMethod
	Mean <- x$Mean
	Nobs <- x$Nobs
	Nrm  <- x$Nrm
	balanced <- x$balanced
	NegVCmsg <- ifelse(is.null(x$NegVCmsg), "", x$NegVCmsg)
	
	if(x$Type == "Linear Model")
	{
		if(!"skipHeading" %in% names(list(...)))
			cat("\n\nAnalysis of Variance Table:\n---------------------------\n\n")
		tmp <- x$aov.org
		tmp$DF <- round(tmp$DF)
		printCoefmat(tmp, digits=digits, dig.tst=digits, has.Pvalue=TRUE, P.values=5, 
				na.print="", zap.ind=3, tst.ind=4, cs.ind=NULL)
		cat("\n")
		
	}
	else
	{
		if(x$EstMethod == "ANOVA")
		{
			MM <- x$Type == "Mixed Model" 
			if(!"skipHeading" %in% names(list(...)) && !MM)
				cat("\n\nResult Variance Component Analysis:\n-----------------------------------\n\n")
			if(!"skipHeading" %in% names(list(...))&& MM)
			{
				cat("\n\nANOVA-Type Estimation of Mixed Model:\n--------------------------------------\n\n")
			}
			if(MM)
			{
				fe <- fixef(x)[,"Estimate", drop=FALSE]
				nam <- rownames(fe)
				fe <- c(fe)
				names(fe) <- nam
			}
			rn <- rownames(x$aov.tab)
			cn <- colnames(x$aov.tab)
			
			#		if(!MM)
			NegVCind <- which(x$VCoriginal < 0) + 1                    	# "VCoriginal" does not contain total variance --> +1
			mat <- matrix(x$aov.tab, nrow=length(rn), dimnames=list(rn, cn))
			mat <- apply(mat, 1:2, round, digits=digits)
			mat <- cbind(Name=rn, mat)
			if("CV" %in% colnames(mat))
				colnames(mat)[which(colnames(mat) == "CV")] <- "CV[%]"
			rownames(mat) <- 1:nrow(mat)
			if(NegVCmsg != "")
			{
				mat[NegVCind, c("VC", "%Total", "SD", "CV[%]")] <- paste(mat[NegVCind, "VC"], "*", sep="")
			}
			mat <- apply(mat, 1:2, function(x) ifelse(x==NaN, NA, x))
			if(MM)
			{
				cat("\t[Fixed Effects]\n\n")
				print(round(fe, digits))
				cat("\n\n\t[Variance Components]\n\n")
			}
			print(mat, na.print="", quote=FALSE, digits=digits)
		}
		else
		{
			MM <- x$Type == "Mixed Model" 
			if(!"skipHeading" %in% names(list(...)) && !MM)
				cat("\n\nResult Variance Component Analysis:\n-----------------------------------\n\n")
			if(!"skipHeading" %in% names(list(...))&& MM)
			{
				cat("\n\nREML-Estimation of Mixed Model:\n-------------------------------\n\n")
			}
			if(MM)
			{
				fe <- fixef(x)[,"Estimate", drop=FALSE]
				nam <- rownames(fe)
				fe <- c(fe)
				names(fe) <- nam
			}
			rn <- rownames(x$aov.tab)
			cn <- colnames(x$aov.tab)
			mat <- as.matrix(x$aov.tab)
			mat <- apply(mat, 1:2, round, digits=digits)
			mat <- cbind(Name=rn, mat)
			rownames(mat) <- 1:nrow(mat)
			
			if(MM)
			{
				cat("\t[Fixed Effects]\n\n")
				print(round(fe, digits))
				cat("\n\n\t[Variance Components]\n\n")
			}
			
			print(mat, na.print="", quote=FALSE, digits=digits)
		}
	}
	
	cat("\nMean:", round(Mean, digits), 
			paste("(N = ",Nobs, ifelse(is.null(Nrm), "", paste(", ", Nrm, " observations removed due to missing data", sep="")), ")", sep=""), 
			"\n\nExperimental Design:", balanced," |  Method:", x$EstMethod)
	
	if(NegVCmsg != "" && x$Type != "Linear Model")                # there were VCs set to zero
		cat(" |", NegVCmsg, "| adapted MS used for total DF")
	cat("\n\n")
}


#' Inferential Statistics for VCA-Results.
#' 
#' Function \code{VCAinference} constructs one- and two-sided confidence intervals, and performs Chi-Squared tests for total
#' and error variance against claimed values for 'VCA' objects.
#' 
#' This function computes confidence intervals (CI) for variance components (VC), standard deviations (SD)
#' and coefficients of variation (CV). VCs 'total' and 'error' can be tested against claimed values specifying parameters
#' 'total.claim' and 'error.claim'. One can also specify claim-values in terms of SD or CV (see \code{claim.type}).\cr
#' Confidence intervals for VCs are constructed either following the same rules as in SAS 9.2 PROC MIXED with option 'method=type1'
#' (ci.method="sas") or using Satterthwaite methodology throughout (ci.method="satterthwaite"). In the former approach
#' for VC total and error, which are constrained to be \eqn{>= 0}, CIs are based on the Chi-Squared distribution. Degrees of freedom
#' (DF) for total variance are approximated using the Satterthwaite approximation (which is not available in either SAS procedure).
#' For all other VCs, the CI is \eqn{[sigma^2-QNorm(alpha/2)*SE(sigma^2); sigma^2+QNorm(1-alpha/2)*SE(sigma^2)]}, where QNorm(x) indicates the x-quantile of 
#' the standard normal distribution. The second method approximates DFs for all VCs using the Satterthwaite approximation and CIs are
#' based on the corresponding Chi-Squared distribution for all VCs (see examples). 
#' Note that in the computation of the covariance-matrix of the VCs, the estimated VCs will be used. If these are requested to be set to 0 
#' (\code{NegVC=FALSE} in \code{\link{anovaVCA}}), the result might not be conformable with theory given in the first reference. 
#' The validity of this implementation was checked against SAS 9.2 PROC MIXED (method=type1), where VCs are not constrained to be >= 0. 
#' The sampling variances for VCs are obtained assuming normality throughout based on \eqn{Var(\sigma^{2} = C^{-1} * Var(m_{SS} * (C^{-1})^{T}))}{Var(sigma^2) = Ci * Var(SS) * Ci'}, 
#' where \eqn{C^{-1}}{Ci} is the inverse of the coefficient matrix equating observed Sum of Squares (SS)
#' to their expected values, and \eqn{(C^{-1})^{T}}{Ci'} indicating the transpose of \eqn{C^{-1}}{Ci} (see Searle et al. 1992, pg. 176).
#' 
#' An input \code{VCA}-object can be in one of three states:\cr 
#' \cr
#'   State (1) corresponds to the situation, where all VC > 0.\cr 
#'   State (2) corresponds to the situation, where at least one VC < 0.\cr
#'   State (3) corresponds to situations, where negative VC estimates occured but were set to 0, i.e. \code{NegVC=FALSE} - the Default.\cr
#' 
#' State (2) occurs when parameter \code{NegVC} was set to TRUE in \code{\link{anovaVCA}}, state (3) represents the default-setting in 
#' function \code{\link{anovaVCA}}. If a \code{VCA}-object is in state (1), parameter \code{excludeNeg} has no effect (there are no negative VCs),
#' only parameter \code{constrainCI} is evaluated. For \code{VCA}-objects in state(2), \code{constrainCI} has no effect, because constraining
#' CIs for unconstrained VCs makes no sense. State (3) forces parameter \code{constrainCI} to be set to TRUE and one can only choose whether to
#' exclude CIs of negative VC estimates or not. Whenever VCs have to be constrained, it is straight forward to apply constraining also to any
#' CI. Note that situations outlined above only occur when parameter \code{VarVC} is set to TRUE, which causes estimation of the covariance-matrix
#' of variance components. The default is only to compute and report CIs for total and error variance, which cannot become negative.
#' 
#' 
#' @param obj               (object) of class 'VCA' or, alternatively, a list of 'VCA' objects, where all other argument can be 
#' 							specified as vectors, where the i-th vector element applies to the i-th element of 'obj' (see examples) 
#' @param alpha             (numeric) value specifying the significance level for \eqn{100*(1-alpha)}\% confidence intervals.
#' @param total.claim       (numeric) value specifying the claim-value for the Chi-Squared test for the total variance (SD or CV, see \code{claim.type}).
#' @param error.claim       (numeric) value specifying the claim-value for the Chi-Squared test for the error variance (SD or CV, see \code{claim.type}).
#' @param claim.type        (character) one of "VC", "SD", "CV" specifying how claim-values have to be interpreted:\cr
#'                          "VC" (Default) = claim-value(s) specified in terms of variance(s),\cr
#'                          "SD" = claim-values specified in terms of standard deviations (SD),\cr
#'                          "CV" = claim-values specified in terms of coefficient(s) of variation (CV)
#'                          and are specified as percentages.\cr
#'                          If set to "SD" or "CV", claim-values will be converted to variances before applying the Chi-Squared test (see examples).
#' @param VarVC             (logical) TRUE = if element "Matrices" exists (see \code{\link{anovaVCA}}), the covariance
#'                          matrix of the estimated VCs will be computed (see \code{\link{vcovVC}}, which is used in CIs for 
#' 							intermediate VCs if 'method.ci="sas"'. 
#'                          Note, this might take very long for larger datasets, since there are many matrix operations involved. 
#'                          FALSE (Default) = computing covariance matrix of VCs is omitted, as well as CIs for intermediate VCs.
#' @param excludeNeg        (logical) TRUE = confidence intervals of negative variance estimates will not be reported. \cr
#'                          FALSE = confidence intervals for all VCs will be reported including those with negative VCs.\cr
#'                          See the details section for a thorough explanation.
#' @param constrainCI       (logical) TRUE = CI-limits for all variance components are constrained to be >= 0.\cr
#'                          FALSE = unconstrained CIs with potentially negative CI-limits will be reported.\cr
#'                          which will preserve the original width of CIs.
#'                          See the details section for a thorough explanation.
#' @param ci.method			(character) string or abbreviation specifying which approach to use for computing confidence intervals of variance components (VC).
#' 							"sas" (default) uses Chi-Squared based CIs for total and error and normal approximation for all other VCs (Wald-limits, option "NOBOUND"
#' 							in SAS PROC MIXED); "satterthwaite" will approximate DFs for each VC using the Satterthwaite approach (see \code{\link{SattDF}} for models
#' 							fitted by ANOVA) and all Cis are based on the Chi-Squared distribution. This approach is conservative but avoids negative values for the lower bounds. 
#' @param quiet				(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' @note                    Original CIs will always be available independent of parameter-settings of \code{excludeNeg} and
#'                          \code{constrainCI}. Original CIs are stored in attribute "CIoriginal" of the returned 'VCAinference'-object, e.g.
#' 							'attr(obj$ConfInt$SD$OneSided, "CIoriginal")' or 'attr(obj$ConfInt$CV$TwoSided, "CIoriginal")'.
#' 
#' @return  (VCAinference) object, a list with elements:
#'          \item{ChiSqTest}{(data.frame) with results of the Chi-Squared test}
#'          \item{ConfInt}{(list) with elements \code{VC}, \code{SD}, \code{CV}, all lists themselves containing (data.frame) objects \code{OneSided} and \code{TwoSided}}
#'          \item{VCAobj}{(VCA) object specified as input, if \code{VarVC=TRUE}, the 'aov.tab' element will have an extra column "Var(VC)" storing variances of VC-estimates"}
#'     
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @references 
#' 
#' Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components., Wiley New York
#' 
#' Burdick, R., Graybill, F. (1992), Confidence Intervals on Variance Components. Marcel Dekker, Inc.
#' 
#' Satterthwaite, F.E. (1946), An Approximate Distribution of Estimates of Variance Components., 
#' Biometrics Bulletin 2, 110-114
#' 
#' @seealso \code{\link{print.VCAinference}}, \code{\link{getVCvar}}, \code{\link{anovaVCA}}
#' 
#' @examples 
#' \dontrun{
#' 
#' # load data (CLSI EP05-A2 Within-Lab Precision Experiment) 
#' data(dataEP05A2_1)
#' 
#' # perform (V)variance (C)component (A)nalysis (also comute A-matrices)
#' res <- anovaVCA(y~day/run, dataEP05A2_1)
#' 
#' # get confidence intervals for total and error (VC, SD, CV)
#' VCAinference(res)
#' 
#' # additionally request CIs for all other VCs; default is to constrain 
#' # CI-limits to be >= 0
#' # first solve MME
#' res <- solveMME(res)
#' VCAinference(res, VarVC=TRUE)
#' 
#' # now using Satterthwaite methodology for CIs
#' VCAinference(res, VarVC=TRUE, ci.method="satt")
#' 
#' # request unconstrained CIs
#' VCAinference(res, VarVC=TRUE, constrainCI=FALSE)
#' 
#' # additionally request Chi-Squared Tests of total and error, default 
#' # is that claim values are specified as variances (claim.type="VC")
#' VCAinference(res, total.claim=4.5, error.claim=3.5)
#' 
#' # perform Chi-Squared Tests, where claim-values are given as SD, 
#' # compare p-values to former example
#' VCAinference(res, total.claim=sqrt(4.5), error.claim=sqrt(3.5), claim.type="SD")
#' 
#' # now using Satterthwaite methodology for CIs
#' VCAinference(res, total.claim=sqrt(4.5), error.claim=sqrt(3.5), 
#' 				claim.type="SD", ci.method="satt")
#' 
#' # now add random error to example data forcing the ANOVA-estimate of the 
#' # day-variance to be negative
#' set.seed(121)
#' tmpData <- dataEP05A2_1
#' tmpData$y <- tmpData$y + rnorm(80,,3)
#' res2 <- anovaVCA(y~day/run, tmpData)
#' 
#' # call 'VCAinference' with default settings
#' VCAinference(res2)
#' 
#' # extract components of the returned 'VCAinference' object
#' inf <- VCAinference(res2, total.claim=12)
#' inf$ConfInt$VC$OneSided			# one-sided CIs for variance components
#' inf$ConfInt$VC$TwoSided			# two-sided CI for variance components
#' inf$ChiSqTest
#' 
#' # request CIs for all VCs, default is to exclude CIs of negative VCs (excludeNeg=TRUE) 
#' # solve MMEs first (or set MME=TRUE when calling anovaVCA)
#' res2 <- solveMME(res2)
#' VCAinference(res2, VarVC=TRUE)
#' 
#' # request CIs for all VCs, including those for negative VCs, note that all CI-limits 
#' # are constrained to be >= 0
#' VCAinference(res2, VarVC=TRUE, excludeNeg=FALSE)
#' 
#' # request unconstrained CIs for all VCs, including those for negative VCS
#' # one has to re-fit the model allowing the VCs to be negative
#' res3 <- anovaVCA(y~day/run, tmpData, NegVC=TRUE, MME=TRUE)
#' VCAinference(res3, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE)
#'  
#' ### use the numerical example from the CLSI EP05-A2 guideline (p.25)
#' data(Glucose)
#' res.ex <- anovaVCA(result~day/run, Glucose)
#' 
#' ### also perform Chi-Squared tests
#' ### Note: in guideline claimed SD-values are used, here, claimed variances are used
#' VCAinference(res.ex, total.claim=3.4^2, error.claim=2.5^2)
#' 
#' 
#' # load another example dataset and extract the "sample_1" subset
#' data(VCAdata1)
#' sample1 <- VCAdata1[which(VCAdata1$sample==1),]
#' 
#' # generate an additional factor variable and random errors according to its levels
#' sample1$device <- gl(3,28,252)                                      
#' set.seed(505)
#' sample1$y <- sample1$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     
#' 
#' # fit a crossed-nested design with main factors 'lot' and 'device' 
#' # and nested factors 'day' and 'run' nested below, also request A-matrices 
#' res1 <- anovaVCA(y~(lot+device)/day/run, sample1) 
#' 
#' # get confidence intervals, covariance-matrix of VCs, ..., 
#' # explicitly request the covariance-matrix of variance components
#' # solve MMEs first
#' res1 <- solveMME(res1)
#' inf1 <- VCAinference(res1, VarVC=TRUE, constrainCI=FALSE)
#' inf1
#' 
#' # print numerical values with more digits
#' print(inf1, digit=12)
#' 
#' # print only parts of the 'VCAinference' object (see \code{\link{print.VCAinference}})
#' print(inf1, digit=12, what=c("VCA", "VC"))
#' 
#' # extract complete covariance matrix of variance components 
#' # (main diagonal is part of standard output -> "Var(VC"))
#' VarCovVC <- vcovVC(inf1$VCAobj)
#' round(VarCovVC, 12)
#' 
#' # use by-processing and specific argument-values for each level of the by-variable
#' data(VCAdata1)
#' fit.all <- anovaVCA(y~(device+lot)/day/run, VCAdata1, by="sample", NegVC=TRUE)
#' inf.all <- VCAinference(fit.all, total.claim=c(.1,.75,.8,1,.5,.5,2.5,20,.1,1))
#' print.VCAinference(inf.all, what="VC")
#' }

VCAinference <- function(obj, alpha=.05, total.claim=NA, error.claim=NA, claim.type="VC", 
		VarVC=FALSE, excludeNeg=TRUE, constrainCI=TRUE, ci.method="sas",
		quiet=FALSE)
{
	Call <- match.call()
	
	if(is.list(obj) && class(obj) != "VCA")
	{
		if(!all(sapply(obj, class) == "VCA"))
			stop("Only lists of 'VCA' object are accepted!")
		
		obj.len <- length(obj)
		
		if(!"msgEnv" %in% ls(.GlobalEnv))
			msgEnv <<- new.env(parent=emptyenv())
		
		assign("VCAinference.obj.is.list", TRUE, envir=msgEnv)			# indicate that a list-type object was passed intially
		
		res <- mapply(	FUN=VCAinference, obj=obj, alpha=alpha, 
				total.claim=total.claim, error.claim=error.claim,
				claim.type=claim.type, VarVC=VarVC, excludeNeg=excludeNeg,
				constrainCI=constrainCI, ci.method=ci.method, 
				SIMPLIFY=FALSE)
		names(res) <- names(obj)
		
		if(obj.len == 1)			# mapply returns a list of length 2 in case that length(obj) was equal to 1
			res <- res[1]
		
		rm("VCAinference.obj.is.list", envir=msgEnv)
		
		return(res)
	}	
	
	stopifnot(class(obj) == "VCA")
	ci.method <- match.arg(ci.method, c("sas", "satterthwaite"))
	
	MM <- obj$Type == "Mixed Model"					# only exists for mixed models
	if(MM)
		VCs   <- obj$Matrices$VCall
	else
		VCs <- obj$aov.tab[-1,"VC"]		 			# for random models directly from ANOVA-table
	
	claim.type <- match.arg(claim.type, c("VC", "SD", "CV"))
	
	Mean <- obj$Mean
	Nvc  <- nrow(obj)
	EstMethod <- obj$EstMethod
	
	if( all(obj$aov.tab[,"VC"] > 0) )
		VCstate <- 1                               
	else if( any(obj$aov.tab[,"VC"] < 0) )
		VCstate <- 2
	else
	{
		if(obj$NegVCmsg != "")
			VCstate <- 3                            # negative VCs set to zero                               
		else                                                                
			VCstate <- 1                            # zero-VC estimated as such                                      
	}
	
	if( !is.null(obj$Matrices) && VarVC )
	{
		if(is.null(obj$VarCov))						# solve mixed model equations first
		{
			obj  <- solveMME(obj)
			
			Lmat  <- obj$Matrices                               # different matrices needed for VCA
			
			VCvar <- vcovVC(obj, method=obj$VarVC.method) 		# get variance-covariance matrix of VCs (p.176); do not pass total VC
			
			NegVCmsg    <- obj$NegVCmsg
			VCoriginal  <- obj$VCoriginal
			Nobs        <- obj$Nobs
			Nrm         <- obj$Nrm
			balanced    <- obj$balanced
			
			obj$aov.tab <- cbind(obj$aov.tab, "Var(VC)"=c(NA, diag(VCvar)))  
			
			class(obj) <- "VCA"
			obj$NegVCmsg   <- NegVCmsg
			obj$VCoriginal <- VCoriginal
			obj$VarCov     <- VCvar                                    # store variance-covariance matrix of variance components
			obj$Mean       <- Mean
			obj$Nobs       <- Nobs
			obj$Nrm        <- Nrm
			obj$balanced   <- balanced
			
			nam0 <- deparse(Call$obj)
			nam1 <- sub("\\[.*", "", nam0)
			
			if(length(nam1) == 1 && nam1 %in% names(as.list(.GlobalEnv)))		# obj is not function call
			{
				expr <- paste(nam0, "<<- obj")						# update object missing MME results
				eval(parse(text=expr))	
			}
			else	# warning only if not called on list of VCA-objects
			{
				if( !"VCAinference.obj.is.list" %in% names(as.list(msgEnv)) && !quiet )
					warning("Mixed model equations solved locally. Results could not be assigned to object!")
			}
		}
	}
	
	if(!is.na(total.claim) && nrow(obj$aov.tab) == 1)                               # if error is the only VC no total variance exists (or is equal)
		total.claim <- NA
	
	if(nrow(obj$aov.tab) == 1)
		Nvc <- 1
	else
		Cind <- which(rownames(obj$aov.tab) %in% c("total", "error"))
	
	if(!is.na(total.claim))
	{
		if(is.numeric(total.claim))
		{
			if( (is.na(total.claim) || total.claim <= 0 ) && !quiet)
				warning("Parameter 'total.claim' is not correctly specified! Chi-Squared test for total precision is omitted!")
		}
		else
		{
			if(!quiet)
				warning("Parameter 'total.claim' is not correctly specified! Chi-Squared test for total precision is omitted!")
		}
	}
	if(!is.na(error.claim))
	{
		if(is.numeric(error.claim))
		{
			if( (is.na(error.claim) || error.claim <= 0) && !quiet)
				warning("Parameter 'error.claim' is not correctly specified! Chi-Squared test for error precision (repeatability) is omitted!")
		}
		else
		{
			if(!quiet)
				warning("Parameter 'error.claim' is not correctly specified! Chi-Squared test for error precision (repeatability) is omitted!")
		}
	}
	
	nam <- rownames(obj$aov.tab)
	Nvc <- length(nam)                                  # number of variance components including total variance
	
	# Chi-Squared tests
	
	CStest <- data.frame(Name=nam, "Claim"=rep(NA,Nvc))
	CStest$"ChiSq value" <- rep(NA, Nvc)
	CStest$"Pr (>ChiSq)" <- rep(NA, Nvc)
	
	if(!is.na(total.claim))
	{
		if(claim.type == "SD")
			total.claim <- total.claim^2                # now the Chi-Squared Test for variances can be used
		
		if(claim.type =="CV")
			total.claim <- (total.claim * Mean/100)^2   # re-caluculate variance from CV
		
		CStest$"Claim"[1] <- total.claim
		
		CStest$"ChiSq value"[1] <- obj$aov.tab[1, "VC"]*obj$aov.tab[1, "DF"]/total.claim 
		
		CStest$"Pr (>ChiSq)"[1] <- pchisq(q=CStest$"ChiSq value"[1], df=obj$aov.tab[1, "DF"], lower.tail=TRUE)
		
		if(claim.type == "SD")
			CStest$"Claim"[1] <- sqrt(total.claim)
		
		if(claim.type == "CV")
			CStest$"Claim"[1] <- sqrt(total.claim)*100/Mean
	}
	
	if(!is.na(error.claim))
	{
		if(claim.type == "SD")
			error.claim <- error.claim^2                # now the Chi-Squared Test for variances can be used
		
		if(claim.type =="CV")
			error.claim <- (error.claim * Mean/100)^2   # re-caluculate variance from CV
		
		CStest$"Claim"[Nvc]       <- error.claim
		CStest$"ChiSq value"[Nvc] <- obj$aov.tab[Nvc, "VC"]*obj$aov.tab[Nvc, "DF"]/error.claim 
		CStest$"Pr (>ChiSq)"[Nvc] <- pchisq(q=CStest$"ChiSq value"[Nvc], df=obj$aov.tab[Nvc, "DF"], lower.tail=TRUE)
		
		if(claim.type == "SD")
			CStest$"Claim"[Nvc] <- sqrt(error.claim)
		
		if(claim.type == "CV")
			CStest$"Claim"[Nvc] <- sqrt(error.claim)*100/Mean
	}
	
	rownames(CStest) <- CStest$Name
	
	# CIs on VCs, SDs and CVs
	
	indCS <- which(rownames(obj$aov.tab) %in% c("total", "error"))				# Chi-Squred dist VCs
	indN  <- which(!rownames(obj$aov.tab) %in% c("total", "error"))				# Normal dist VCs
	
	#####################
	### two-sided CIs ###
	
	CI_VC <- data.frame(Name=nam)
	CI_VC$LCL <- numeric(nrow(obj$aov.tab))
	CI_VC$UCL <- numeric(nrow(obj$aov.tab))
	
	# Chi-Squared CIs
	
	lower.qchisq <- qchisq(p=1-alpha/2, df=obj$aov.tab[indCS, "DF"])
	upper.qchisq <- qchisq(p=alpha/2, df=obj$aov.tab[indCS, "DF"])
	CI_VC$LCL[indCS] <- obj$aov.tab[indCS, "VC"]*obj$aov.tab[indCS, "DF"]/lower.qchisq
	CI_VC$UCL[indCS] <- obj$aov.tab[indCS, "VC"]*obj$aov.tab[indCS, "DF"]/upper.qchisq
	
	# CIs based on Normal-Distribution (ci.method="sas") or Chi-Squared (ci.method="satterthwaite")
	
	if(ci.method == "sas" && "Var(VC)" %in% colnames(obj$aov.tab))		# variance-covariance matrix of VCs required
	{
		CI_VC$LCL[indN] <- obj$aov.tab[indN, "VC"]+qnorm(alpha/2)*sqrt(obj$aov.tab[indN, "Var(VC)"])
		CI_VC$UCL[indN] <- obj$aov.tab[indN, "VC"]+qnorm(1-alpha/2)*sqrt(obj$aov.tab[indN, "Var(VC)"])
		DFs <- NULL
	}
	else if(ci.method == "satterthwaite")
	{
		if(obj$EstMethod == "ANOVA")
		{
			rind	<- obj$Matrices$rf.ind						# indices of random terms in the original ANOVA-table
			aov.tab <- if(obj$Type == "Random Model")			# original ANOVA table
						obj$aov.tab[-1,]					# remove row for total 
					else	
						obj$aov.org
			MS  	<- aov.tab[rind, "MS"]
			Ci  	<- getMat(obj, "Ci.MS")[rind, rind]
			DF  	<- aov.tab[rind, "DF"]
			sDF 	<- SattDF(MS, Ci, DF, type="individual") 	# satterthwaite DFs
			ind 	<- 1:(length(sDF)-1)						# without row for error 
			DFs 	<- c(obj$aov.tab[1,"DF"], sDF)				# total DF
		}
		else													# fitted by REML, only Satterthwaite DF exist
		{
			DFs <- obj$aov.tab[,"DF"]
			sDF <- DFs[-1]
			ind <- 1:(length(sDF)-1)
		}
		
		CI_VC$DF <- DFs										# add Satterthwaite DFs
		CI_VC <- CI_VC[,c(1,4,2,3)]
		lower.qchisq <- qchisq(p=1-alpha/2, df=sDF[ind])
		upper.qchisq <- qchisq(p=alpha/2, df=sDF[ind])
		CI_VC$LCL[ind+1] <- obj$aov.tab[ind+1, "VC"]*sDF[ind]/lower.qchisq		# lower limit all but total and error
		CI_VC$UCL[ind+1] <- obj$aov.tab[ind+1, "VC"]*sDF[ind]/upper.qchisq		# upper limit  - " -
	}
	else
	{
		CI_VC$LCL[indN] <- NA
		CI_VC$UCL[indN] <- NA
		DFs <- NULL
	}
	
	rownames(CI_VC) <- CI_VC$Name
	attr(CI_VC, "CIoriginal") <- CI_VC                      # keep original CI limits as attribute
	
	if(VCstate == 1)                                        # excludeNeg not evaluated, since all VCs positive
	{
		if( constrainCI && any(CI_VC$LCL < 0, na.rm=TRUE))  # any negative CI-limits?
		{
			LCLind <- which(CI_VC$LCL < 0)
			UCLind <- which(CI_VC$UCL < 0)                  # keep info about constrained LCLs
			
			attr(CI_VC, "LCLconstrained") <- LCLind
			
			CI_VC$LCL[LCLind] <- 0                          # set negative LCL to 0
			
			if(length(UCLind) > 0)                          # set negative UCL to 0 
			{
				CI_VC$UCL[UCLind] <- 0
				attr(CI_VC, "UCLconstrained") <- UCLind     # keep info about constrained UCLs
			}
		}
	}
	if(VCstate == 2)                                        # constrainCI not evaluated since negative VCs explicitly allowed
	{
		if( excludeNeg )                                    # CIs for negative VC are excluded
		{    
			NegInd <- which(obj$aov.tab[,"VC"] < 0)
			CI_VC$LCL[NegInd] <- NA                         # there has to be at least one LCL < 0
			CI_VC$UCL[NegInd] <- NA
		}
	}
	if(VCstate == 3)
	{
		if( excludeNeg )
		{
			NegInd <- which(obj$VCoriginal < 0)+1  			# which VCs were set to 0, i.e. were originally < 0
			CI_VC$LCL[NegInd] <- NA                         # there has to be at least one LCL < 0
			CI_VC$UCL[NegInd] <- NA
		}
		
		if(any(CI_VC$LCL < 0, na.rm=TRUE))                  # constrainCI implicitly set to 0, since negative VCs set to zero automatically results in constraining CIs
		{
			LCLind <- which(CI_VC$LCL < 0)
			UCLind <- which(CI_VC$UCL < 0)                  # keep info about constrained LCLs
			
			attr(CI_VC, "LCLconstrained") <- LCLind
			
			CI_VC$LCL[LCLind] <- 0                          # set negative LCL to 0
			
			if(length(UCLind) > 0)                          # set negative UCL to 0 
			{
				CI_VC$UCL[UCLind] <- 0
				attr(CI_VC, "UCLconstrained") <- UCLind     # keep info about constrained UCLs
			}
		}
	}
	
	# test whether point-estimates within CI or not
	CIwarning <- FALSE
	CIout <- NULL
	CIind <- which(obj$aov.tab[,"VC"] < CI_VC$LCL | obj$aov.tab[,"VC"] > CI_VC$UCL)
	
	if(length(CIind) > 0)
	{
		CI_VC$LCL[CIind] <- CI_VC$UCL[CIind] <- NA									# set them NA
		CIwarning <- TRUE
		CIout <- c(CIout, rownames(obj$aov.tab)[CIind])
	}
	
	
	# SD and CV computed from VC      NOTE: sqrt o.k.??????????
	
	CI_SD <- data.frame(Name=nam)
	CI_SD$DF <- DFs							# add Satterthwaite DFs if available
	Sign <- sign(CI_VC$LCL)
	CI_SD$LCL <- sqrt(abs(CI_VC$LCL))*Sign
	Sign <- sign(CI_VC$UCL)
	CI_SD$UCL <- sqrt(abs(CI_VC$UCL))*Sign
	rownames(CI_SD) <- CI_SD$Name
	
	if(any(is.na(obj$aov.tab[,"SD"])))
		CI_SD[which(is.na(obj$aov.tab[,"SD"])),2:ncol(CI_SD)] <- NA 
	
	CI_CV <- data.frame(Name=nam)
	CI_CV$DF <- DFs							# add Satterthwaite DFs if available
	CI_CV$LCL <- CI_SD$LCL * 100/Mean
	CI_CV$UCL <- CI_SD$UCL * 100/Mean
	rownames(CI_CV) <- CI_CV$Name
	
	CI_two_sided <- list(CI_VC=CI_VC, CI_SD=CI_SD, CI_CV=CI_CV)
	
	# one-sided CIs
	
	CI_VC <- data.frame(Name=nam)
	CI_VC$DF <- DFs							# add Satterthwaite DFs if available
	CI_VC$LCL <- numeric(nrow(obj$aov.tab))
	CI_VC$UCL <- numeric(nrow(obj$aov.tab))
	
	# CIs based on Chi-Squared
	
	upper.qchisq <- qchisq(p=alpha,   df=obj$aov.tab[indCS, "DF"])
	lower.qchisq <- qchisq(p=1-alpha, df=obj$aov.tab[indCS, "DF"])
	CI_VC$LCL[indCS] <- obj$aov.tab[indCS, "VC"]*obj$aov.tab[indCS, "DF"]/lower.qchisq
	CI_VC$UCL[indCS] <- obj$aov.tab[indCS, "VC"]*obj$aov.tab[indCS, "DF"]/upper.qchisq
	rownames(CI_VC) <- CI_VC$Name
	
	# CIs based on Normal Distribution or Chi-Squared 
	
	if(ci.method == "sas" && "Var(VC)" %in% colnames(obj$aov.tab))	
	{
		CI_VC$UCL[indN] <- obj$aov.tab[indN, "VC"]+qnorm(1-alpha)*sqrt(obj$aov.tab[indN, "Var(VC)"])
		CI_VC$LCL[indN] <- obj$aov.tab[indN, "VC"]+qnorm(  alpha)*sqrt(obj$aov.tab[indN, "Var(VC)"])
	}
	else if(ci.method == "satterthwaite")
	{
		lower.qchisq <- qchisq(p=1-alpha, df=sDF[ind])
		upper.qchisq <- qchisq(p=alpha, df=sDF[ind])
		CI_VC$LCL[ind+1] <- obj$aov.tab[ind+1, "VC"]*sDF[ind]/lower.qchisq
		CI_VC$UCL[ind+1] <- obj$aov.tab[ind+1, "VC"]*sDF[ind]/upper.qchisq
	}
	else
	{
		CI_VC$UCL[indN] <- NA
		CI_VC$LCL[indN] <- NA
	}
	
	attr(CI_VC, "CIoriginal") <- CI_VC                      # keep original CI limits as attribute
	
	if(VCstate == 1)                                        # excludeNeg not evaluated, since all VCs positive
	{
		if( constrainCI && any(CI_VC$LCL < 0, na.rm=TRUE))  # any negative CI-limits?
		{
			LCLind <- which(CI_VC$LCL < 0)
			UCLind <- which(CI_VC$UCL < 0)                  # keep info about constrained LCLs
			CI_VC$LCL[LCLind] <- 0 
			attr(CI_VC, "LCLconstrained") <- LCLind
			
			if(length(UCLind) > 0)                          # set negative UCL to 0 
			{
				CI_VC$UCL[UCLind] <- 0
				attr(CI_VC, "UCLconstrained") <- UCLind     # keep info about constrained UCLs
			}
		}
	}
	if(VCstate == 2)                                        # constrainCI not evaluated since negative VCs explicitly allowed
	{
		if( excludeNeg )                                    # CIs for negative VC are excluded
		{    
			NegInd <- which(obj$aov.tab[,"VC"] < 0)
			
			CI_VC$UCL[NegInd] <- NA                         # there has to be at least one UCL < 0
			CI_VC$LCL[NegInd] <- NA                         # there has to be at least one UCL < 0
		}
	}
	if(VCstate == 3)
	{
		if( excludeNeg )
		{
			NegInd <- which(obj$VCoriginal < 0) + 1  		# which VCs were set to 0, i.e. were originally < 0
			CI_VC$UCL[NegInd] <- NA                         # there has to be at least one UCL < 0
			CI_VC$LCL[NegInd] <- NA
		}
		else                                                # constrainCI implicitly set to TRUE, because some VCs set to zero, unconstrained CIs would not make sense
		{
			if(any(CI_VC$LCL < 0, na.rm=TRUE))
			{
				LCLind <- which(CI_VC$LCL < 0)              # keep info about constrained LCLs
				CI_VC$LCL[LCLind] <- 0                      # set negative LCL to 0
				attr(CI_VC, "LCLconstrained") <- LCLind
			}
			
#            if(any(CI_VC$UCL < 0, na.rm=TRUE))
#            {
#                UCLind <- which(CI_VC$UCL < 0)              # keep info about constrained UCLs
#                
#                attr(CI_VC, "UCLconstrained") <- UCLind
#                
#                CI_VC$UCL[UCLind] <- 0                      # set negative UCL to 0
#            }
		}
	}
	
	# test whether point-estimates within CI or not
	CIindLCL <- which(obj$aov.tab[,"VC"] < CI_VC$LCL)
	CIindUCL <- which(obj$aov.tab[,"VC"] > CI_VC$UCL)
	
	if(length(CIindLCL) > 0)
	{
		CI_VC$LCL[CIindLCL] <- NA								# set LCL to NA
		CIwarning <- TRUE
		CIout <- c(CIout, rownames(obj$aov.tab)[CIindLCL])
	}
	if(length(CIindUCL) > 0)
	{
		CI_VC$UCL[CIindUCL] <- NA										# set LCL to NA
		CIwarning <- TRUE
		CIout <- c(CIout, rownames(obj$aov.tab)[CIindUCL])
	}
	
	if(CIwarning && !quiet)
	{	
		warning("Point estimate(s) of VC(s) ", paste(paste("'", unique(CIout), "'", sep=""), sep=", ")," found outside of confidence interval, CI was set to 'NA'!")
	}
	# SD and CV computed from VC
	
	CI_SD <- data.frame(Name=nam)
	CI_SD$DF <- DFs											# add Satterthwaite DFs if available
	SignLCL  <- sign(CI_VC$LCL)
	CI_SD$LCL <- sqrt(abs(CI_VC$LCL)) * SignLCL
	SignUCL  <- sign(CI_VC$UCL)
	CI_SD$UCL <- sqrt(abs(CI_VC$UCL)) * SignUCL
	rownames(CI_SD) <- CI_SD$Name
	
	if(any(is.na(obj$aov.tab[,"SD"])))
		CI_SD[which(is.na(obj$aov.tab[,"SD"])),2:ncol(CI_SD)] <- NA 
	
	CI_CV <- data.frame(Name=nam)
	CI_CV$DF <- DFs											# add Satterthwaite DFs if available
	CI_CV$LCL <- CI_SD$LCL * 100/Mean
	CI_CV$UCL <- CI_SD$UCL * 100/Mean
	rownames(CI_CV) <- CI_CV$Name
	
	CI_one_sided <- list(CI_VC=CI_VC, CI_SD=CI_SD, CI_CV=CI_CV)
	
	CIs <- list(VC=list(OneSided=CI_one_sided$CI_VC, TwoSided=CI_two_sided$CI_VC),
			SD=list(OneSided=CI_one_sided$CI_SD, TwoSided=CI_two_sided$CI_SD),
			CV=list(OneSided=CI_one_sided$CI_CV, TwoSided=CI_two_sided$CI_CV))
	
	result <- list(ChiSqTest=CStest, ConfInt=CIs, VCAobj=obj, alpha=alpha)
	class(result) <- "VCAinference"
	attr(result, "EstMethod")   <- EstMethod
	attr(result, "excludeNeg")  <- excludeNeg
	attr(result, "constrainCI") <- constrainCI
	attr(result, "claim.type")  <- claim.type
	attr(result, "ci.method")   <- ci.method
	return(result)
}


#' Standard Print Method for Objects of Class 'VCAinference'.
#' 
#' Formats the list-type objects of class 'VCAinference' for a more comprehensive
#' presentation of results, which are easier to grasp. The default is to show the complete
#' object (VCA ANOVA-table, VC-, SD-, and CV-CIs). Using parameter 'what' allows to
#' restrict the printed output to certain parts. Print-function invisibly returns a matrix
#' or a list of matrices, depending on the values of 'what', i.e. it can be used as for
#' packing the inference-information in one or multiple matrix-objects and extracting it/them.
#' 
#' @param x         (VCAinference) object 
#' @param digits    (integer) number of decimal digits.
#' @param what      (character) one of "all", "VC", "SD", "CV", "VCA" specifying which part of the 'VCA'-object is to be printed.
#' @param ...       additional arguments to be passed to or from methods.
#' 
#' @return invisibly returns sub-elements of 'x' specified via 'what'
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @method print VCAinference
#' @S3method print VCAinference
#' 
#' @seealso \code{\link{VCAinference}}, \code{\link{anovaVCA}}
#' 
#' @examples
#' \dontrun{
#' # load data (CLSI EP05-A2 Within-Lab Precision Experiment) 
#' data(dataEP05A2_1)
#' 
#' # perform ANOVA-estimation of variance components for a nested design
#' res <- anovaVCA(y~day/run, Data=dataEP05A2_1)
#' res
#' inf <- VCAinference(res)
#' inf
#' 
#' # show certain parts and extract them invisibly
#' CVmat <- print(inf, what="CV")
#' CVmat
#' 
#' # show numerical values with more digits
#' print(inf, digit=12)
#' }

print.VCAinference <- function(x, digits=4L, what=c("all", "VC", "SD", "CV", "VCA"), ...)
{
	if(is.list(x) && class(x) != "VCAinference")
	{
		if(!all(sapply(x, class) == "VCAinference"))
			stop("Only lists of 'VCAinference' objects can printed!")
		
		nam <- names(x)
		lst <- list()
		
		for(i in 1:length(x))
		{
			print(nam[i])
			lst[[i]] <- print(x[[i]], digits=digits, what=what)
		}
		return()			# leave function now
	}
	
	stopifnot(class(x) == "VCAinference") 
	claim.type <- attr(x, "claim.type")
	VCAobj <- x$VCAobj
	
	ret <- list()
	
	what <- match.arg(tolower(what), tolower(c("all", "VC", "SD", "CV", "VCA")), several.ok=TRUE)
	
#	MM <- !is.null(attr(x$VCAobj, "FixedEffects"))
	
	if(VCAobj$Type == "Mixed Model")
		cat("\n\n\nInference from Mixed Model Fit\n------------------------------\n\n")
	else if(VCAobj$Type == "Linear Model")
		cat("\n\n\nInference from Linear Model Fit:\n--------------------------------\n\n")
	else
		cat("\n\n\nInference from (V)ariance (C)omponent (A)nalysis\n------------------------------------------------\n\n")
	
	if( any( c("all", "vca") %in% what) )
	{
		if(VCAobj$Type == "Linear Model")
			cat("> ANOVA Table:\n--------------\n\n")
		else
			cat("> VCA Result:\n-------------\n\n")
		print(VCAobj, digits=digits, skipHeading=TRUE)
		cat("\n")
	}    
	
	if( any( c("all", "vc") %in% what) )
	{
		cat("> VC:\n-----\n")
		VC <- as.data.frame(as.matrix(VCAobj$aov.tab)[,"VC", drop=FALSE])
		colnames(VC) <- "Estimate"
		
		if(any(!is.na(x$ChiSqTest[,"Claim"])) && claim.type == "VC")
		{
			VC$Claim <- as.numeric(rep(NA, nrow(VC)))
			VC$"ChiSq" <- as.numeric(rep(NA, nrow(VC)))
			VC$"Pr(>ChiSq)" <- as.numeric(rep(NA, nrow(VC)))
			VC[as.character(x$ChiSqTest[,"Name"]), c("Claim", "ChiSq", "Pr(>ChiSq)")] <- x$ChiSqTest[,-1]
		}
		
		VC$DF <- x$ConfInt$VC$TwoSided$DF
		
		VC$"CI LCL" <-  as.numeric(rep(NA, nrow(VC)))
		VC$"CI UCL" <-  as.numeric(rep(NA, nrow(VC)))
		
		VC[as.character(x$ConfInt$VC$TwoSided[,"Name"]), c("CI LCL", "CI UCL")] <- x$ConfInt$VC$TwoSided[,c("LCL", "UCL")]
		
		VC$"One-Sided LCL" <-  as.numeric(rep(NA, nrow(VC)))
		VC$"One-Sided UCL" <-  as.numeric(rep(NA, nrow(VC)))
		
		VC[as.character(x$ConfInt$VC$OneSided[,"Name"]), c("One-Sided LCL", "One-Sided UCL")] <-x$ConfInt$VC$OneSided[,c("LCL", "UCL")]
		
		VC <- as.matrix(VC)
		VC <- apply(VC, 1:2, round, digits)   
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "LCLconstrained")))
			VC[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"] <- paste(VC[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "UCLconstrained")))
			VC[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"] <- paste(VC[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "LCLconstrained")))
			VC[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"] <- paste(VC[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "UCLconstrained")))
			VC[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"] <- paste(VC[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"], "*", sep="")
		
		print(noquote(VC), na.print="")
		ret$VC <- VC
	}
	
	if( any( c("all", "sd") %in% what) )
	{
		cat("\n> SD:\n-----\n")
		SD <- as.data.frame(as.matrix(VCAobj$aov.tab[,"SD", drop=FALSE]))
		colnames(SD) <- "Estimate"
		
		if(any(!is.na(x$ChiSqTest[,"Claim"])) && claim.type == "SD")
		{
			SD$Claim <- as.numeric(rep(NA, nrow(SD)))
			SD$"ChiSq" <- as.numeric(rep(NA, nrow(SD)))
			SD$"Pr(>ChiSq)" <- as.numeric(rep(NA, nrow(SD)))
			SD[as.character(x$ChiSqTest[,"Name"]), c("Claim", "ChiSq", "Pr(>ChiSq)")] <- x$ChiSqTest[,-1]
		}
		
		SD$DF <- x$ConfInt$SD$TwoSided$DF
		
		SD$"CI LCL" <-  as.numeric(rep(NA, nrow(SD)))
		SD$"CI UCL" <-  as.numeric(rep(NA, nrow(SD)))
		
		SD[as.character(x$ConfInt$SD$TwoSided[,"Name"]), c("CI LCL", "CI UCL")] <-x$ConfInt$SD$TwoSided[,c("LCL", "UCL")]
		
		SD$"One-Sided LCL" <-  as.numeric(rep(NA, nrow(SD)))
		SD$"One-Sided UCL" <-  as.numeric(rep(NA, nrow(SD)))
		
		SD[as.character(x$ConfInt$SD$OneSided[,"Name"]), c("One-Sided LCL", "One-Sided UCL")] <-x$ConfInt$SD$OneSided[,c("LCL", "UCL")]
		
		SD <- as.matrix(SD)
		SD <- apply(SD, 1:2, round, digits)  
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "LCLconstrained")))                 # Note: attribute "LCLconstrained" exists only for VC since SD and CV computed from VC
			SD[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"] <- paste(SD[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "UCLconstrained")))                 # see above
			SD[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"] <- paste(SD[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "LCLconstrained")))                 # see above
			SD[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"] <- paste(SD[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "UCLconstrained")))                 # see above
			SD[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"] <- paste(SD[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"], "*", sep="")
		
		print(noquote(SD), na.print="")
		
		ret$SD <- SD
	}
	
	if( any( c("all", "cv") %in% what) )
	{
		cat("\n> CV[%]:\n--------\n")
		CV <- as.data.frame(as.matrix(VCAobj$aov.tab)[,"CV[%]", drop=FALSE])
		colnames(CV) <- "Estimate"
		
		if(any(!is.na(x$ChiSqTest[,"Claim"])) && claim.type == "CV")
		{
			CV$Claim <- as.numeric(rep(NA, nrow(CV)))
			CV$"ChiSq" <- as.numeric(rep(NA, nrow(CV)))
			CV$"Pr(>ChiSq)" <- as.numeric(rep(NA, nrow(CV)))
			CV[as.character(x$ChiSqTest[,"Name"]), c("Claim", "ChiSq", "Pr(>ChiSq)")] <- x$ChiSqTest[,-1]
		}
		
		CV$DF <- x$ConfInt$CV$TwoSided$DF
		
		CV$"CI LCL" <-  as.numeric(rep(NA, nrow(CV)))
		CV$"CI UCL" <-  as.numeric(rep(NA, nrow(CV)))
		
		CV[as.character(x$ConfInt$CV$TwoSided[,"Name"]), c("CI LCL", "CI UCL")] <-x$ConfInt$CV$TwoSided[,c("LCL", "UCL")]
		
		CV$"One-Sided LCL" <-  as.numeric(rep(NA, nrow(CV)))
		CV$"One-Sided UCL" <-  as.numeric(rep(NA, nrow(CV)))
		
		CV[as.character(x$ConfInt$CV$OneSided[,"Name"]), c("One-Sided LCL","One-Sided UCL")] <-x$ConfInt$CV$OneSided[,c("LCL", "UCL")]
		
		CV <- as.matrix(CV)
		CV <- apply(CV, 1:2, round, digits)   
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "LCLconstrained")))                 # Note: attribute "LCLconstrained" exists only for VC since SD and CV computed from VC
			CV[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"] <- paste(CV[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "UCLconstrained")))                 # see above
			CV[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"] <- paste(CV[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "LCLconstrained")))                 # see above
			CV[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"] <- paste(CV[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "UCLconstrained")))                 # see above
			CV[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"] <- paste(CV[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"], "*", sep="")
		
		print(noquote(CV), na.print="")
		
		ret$CV <- CV
	}
	
	ci.method <- attr(x, "ci.method")
	
	if( any(c("all", "vc", "sd", "cv") %in% what) )
		cat( paste("\n\n", 100*(1-x$alpha), "% Confidence Level  ", ifelse(attr(x, "excludeNeg") && (VCAobj$aov.tab[,"VC"] <= 0), "|  CIs for negative VCs excluded  ", ""),
						ifelse(!is.null(attr(x$ConfInt$VC$TwoSided, "LCLconstrained")), "| * CI-limits constrained to be >= 0", ""), 
						ifelse(ci.method=="sas", "\nSAS PROC MIXED method used for computing CIs", "\nSatterthwaite methodology used for computing CIs"), sep=""), "\n\n")
	
	if(length(ret) > 1)                     # extracting parts of the 'VCAinference' object
		invisible(ret)
	else if(length(ret) == 1)
		invisible(ret[[1]])
	
}


#' Check Whether Design Is Balanced Or Not.
#' 
#' This function is for internal use only. Thus, it is not exported.
#' 
#' The approach taken here is to check whether each cell defined by one level of a factor are all equal or
#' not. Here, data is either balanced or unbalanced, there is no concept of "planned unbalancedness" as
#' discussed e.g. in Searle et al. (1992) p.4. The expanded (simplified) formula is divided into main factors
#' and nested factors, where the latter are interaction terms. The N-dimensional contingency table, N being the
#' number of main factors, is checked for all cells containing the same number. If there are differences, the
#' dataset is classified as "unbalanced". All interaction terms are tested individually. Firstly, a single factor 
#' is generated from combining factor levels of the first (n-1) variables in the interaction term. The last variable
#' occuring in the interaction term is then recoded as factor-object with \eqn{M} levels. \eqn{M} is the number of factor
#' levels within each factor level defined by the first \eqn{(n-1)} variables in the interaction term. This is done to 
#' account for the independence within sub-classes emerging from the combination of the first \eqn{(n-1)} variables.
#' 
#' @param form      (formula) object defining the experimental design.
#' @param Data      (data.frame) containing all variables appearing in 'form'.
#' @param na.rm     (logical) TRUE = delete rows where any is NA, FALSE = NAs are not removed, if there are NAs in the
#'                  response variable and all information in independent variables is available, then only the design is checked.
#' 
#' 
#' @return (logical) TRUE if data is balanced, FALSE if data is unbalanced (according to the definition of balance used)
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' 
#' \dontrun{
#' data1 <- data.frame(site=gl(3,8), lot=factor(rep(c(2,3,1,2,3,1), 
#'                     rep(4,6))), day=rep(1:12, rep(2,12)), y=rnorm(24,25,1))
#' 
#' # not all combinations of 'site' and 'lot' in 'data1'
#' 
#' VCA:::isBalanced(y~site+lot+site:lot:day, data1)
#' 
#' # balanced design for this model
#' 
#' VCA:::isBalanced(y~lot+lot:day, data1)
#' 
#' # gets unbalanced if observation is NA
#' 
#' data1[1,"y"] <- NA
#' VCA:::isBalanced(y~lot+lot:day, data1)
#' VCA:::isBalanced(y~lot+lot:day, data1, FALSE)
#' }

isBalanced <- function(form, Data, na.rm=TRUE)
{
	form <- terms(form, simplify=TRUE, keep.order=TRUE)
	if(length(attr(form, "factors")) == 0)
		return(TRUE)
	

	
	rn   <- rownames(attr(form, "factors"))
	
	for(i in 2:length(rn))
	{
		if(!is.numeric(Data[[rn[i]]]))				# for all non-numeric variables
		{
			Data[[rn[i]]] <- factor(as.character(Data[[rn[i]]]))			# get rid of non-existing factor levels (e.g. in data sub-sets)
		}
	}
	
	if(na.rm)
	{
		stopifnot(attr(form, "response") == 1)
		resp <- rn[1]
		Data <- Data[,rn]                           # only used variable appearing in the formula
		Data <- na.omit(Data)
	}
	
	fac  <- attr(form, "term.labels")
	
	mainFac  <- character()
	balanced <- TRUE
	
	for(i in 1:length(fac))
	{
		if(!balanced)
			break
		
		tmp  <- unlist(strsplit(fac[i], ":"))		# split interaction terms into variables
		Nvar <- length(tmp) 
		
		if( Nvar == 1)                              # main factor
		{
			mainFac <- c(mainFac, tmp)       
		} 
		else                                        # nested factors ("a:b" interpreted as b nested in a, since no possible main factor "b" is taken into account)
		{
			tmp.df <- data.frame(var1=apply(Data[,tmp[-Nvar], drop=FALSE], 1, function(x){
								return(paste(paste(letters[1:length(x)], x, sep=""), collapse=""))
							}))
			var2    <- character()
			var1Lev <- unique(tmp.df$var1)
			
			for(j in 1:length(var1Lev))             # start testing each term with nested factors ("a:b")
			{
				ind  <- which(tmp.df$var1 == var1Lev[j])
				var2 <- c(var2, as.factor(as.integer(Data[ind, tmp[Nvar]]))) 
			}
			tmp.df$var2 <- var2
			tmp.tab     <- table(tmp.df)           	# generate contingency table -> ...
			
			balanced <- all(c(tmp.tab) == tmp.tab[1,1])		# ... check for equal numbers in cells of the contingency table
		}
	}
	
	if(length(mainFac) > 0 && balanced)           	# check all main factors if no evidence of unbalncedness
	{
		tmp.tab  <- table(Data[,mainFac])
		tmp.df   <- as.data.frame(tmp.tab)
		balanced <- all(tmp.df[,"Freq"] == tmp.df[1, "Freq"])
	}
	
	return(balanced)
}


#' Compute the Trace of a Matrix.
#' 
#' Function computes the sum of main-diagonal elements of a squared matrix.
#' 
#' @param x     	(matrix, Matrix) object
#' @param quiet		(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' @return (numeric) value, the trace of the matrix

Trace <- function(x, quiet=FALSE)
{
	if(ncol(x) != nrow(x) && !quiet)
		warning("Matrix is not quadratic!")
	return(sum(diag(as.matrix(x))))
}



#' Standard 'as.matrix' Method for 'VCA' S3-Objects.
#' 
#' @param x         (VCA) object 
#' @param ...       additional arguments to be passed to or from methods.
#' 
#' @return (matrix) equal to x$aov.tab with additional attributes "Mean" and "Nobs"
#' 
#' @method as.matrix VCA
#' @S3method as.matrix VCA
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' fit <- anovaVCA(y~day/run, dataEP05A2_1)
#' as.matrix(fit)
#' }

as.matrix.VCA <- function(x, ...)
{
	Mean <- x$Mean
	Nobs <- x$Nobs
	rn <- rownames(x$aov.tab)
	cn <- colnames(x$aov.tab)
	mat <- matrix(x$aov.tab, nrow=length(rn), ncol=length(cn), dimnames=list(rn, cn))
	attr(mat, "Mean") <- Mean
	attr(mat, "Nobs") <- Nobs
	return(mat)
}



#' ANOVA-Type Estimation of Variance Components for Random Models.
#' 
#' This function equates observed ANOVA Type-I sums of squares (\eqn{SS}) to their expected values and solves the resulting system of linear equations
#' for variance components. 
#' 
#' For diagnostics, a key parameter is "precision", i.e. the accuracy of a quantification method influenced by varying sources of random error. 
#' This type of experiments is requested by regulatory authorities to proof the quality of diagnostic tests, e.g. quantifying intermediate
#' precision according to CLSI guideline EP5-A2/A3. No, fixed effects are allowed besides the intercept. 
#' Whenever fixed effects are part of the model to be analyzed, use function \code{\link{anovaMM}} instead.
#' 
#' Function \code{anovaVCA} is tailored for performing Variance Component Analyses (VCA) for random models, assuming all VCs as factor variables, i.e. their levels
#' correspond to distinct columns in the design matrix (dummy variables). Any predictor variables are automatically converted to factor variables, since continuous
#' variables may not be used on the right side of the formula 'form'.
#' 
#' ANOVA \eqn{SS} are either computed employing the SWEEP-operator (Goodnight 1979, default) or by finding matrices \eqn{A_i}{A_i} expressing them as quadratic forms 
#' in \eqn{y} as \eqn{ss_i = y^{T}A_{i}y}{ss_i = y' * A_i * y}. Matrices \eqn{A_i} are also used to compute the variance-coveriance matrix of variance components (VC)
#' according to Searle et al. (1992) which corresponds to \code{VarVC.method="scm"}. Whenever the SWEEP-operator is used, which is way faster and therefore the default method,
#' the approximation according to Giesbrecht and Burns (1985) is automatically set (\code{VarVC.method="gb"}). 
#' 
#' Function \code{anovaVCA} represents a special form of the "method of moments" approach applicable to arbitrary random models either balanced or unbalanced.
#' The system of linear equations, which is built from the ANOVA Type-I sums of squares, is closely related to the method used 
#' by SAS PROC VARCOMP, where ANOVA mean squares (\eqn{MS}) are used (see \code{\link{getCmatrix}}). The former can be written as \eqn{ss = C * s}
#' and the latter as \eqn{ms = D * s}, where \eqn{C} and \eqn{D} denote the respective coefficient matrices, \eqn{s} the column-vector
#' of variance components (VC) to be estimated/predicted, and \eqn{ss} and \eqn{ms} the column vector of ANOVA sum of squares, respectively, mean squares. 
#' Mutliplying element \eqn{d_ij} of matrix \eqn{D} by element \eqn{c_in} of matrix \eqn{C} (\eqn{i,j = 1,...,n}), results in 
#' matrix \eqn{C}. Thus, \eqn{C} can easily be converted to \eqn{D} by the inverse operation. Matrix \eqn{D} is used to estimate
#' total degrees of freedom (DF) according to Satterthwaite (1946).
#' 
#' Both methods for computing ANOVA Type-I \eqn{SS} are much faster than fitting the linear model via \code{\link{lm}} and calling function \code{\link{anova}} on the 'lm' object
#' for complex models, where complex refers to the number of columns of the design matrix and the degree of unbalancedness. Degrees of freedom for the \eqn{i}-th term are obtained 
#' by function \code{\link{anovaDF}} in case of \code{VarVC.method="scm"}. Otherwise, \eqn{DF} are directly derived from the SWEEP-operator as the number of linearly independent
#' columns of the partial design matrix corresponding to a specific \eqn{VC}.
#' 
#' @param form          (formula) specifying the model to be fit, a response variable left of the '~' is mandatory
#' @param Data          (data.frame) storing all variables referenced in 'form'
#' @param by			(factor, character) variable specifying groups for which the analysis should be performed individually,
#' 						i.e. by-processing
#' @param NegVC         (logical) FALSE = negative variance component estimates (VC) will be set to 0 and they will not contribute to the total variance 
#'                      (as done in SAS PROC NESTED, conservative estimate of total variance). The original ANOVA estimates can be found in element 'VCoriginal'. 
#'                      The degrees of freedom of the total variance are based on adapted mean squares (MS), i.e. adapted MS are computed as \eqn{D * VC}, where VC is 
#'                      the column vector with negative VCs set to 0. \cr
#' 						TRUE = negative variance component estimates will not be set to 0 and they will contribute to the total variance (original definition of the total variance).
#' @param SSQ.method	(character) string specifying the method used for computing ANOVA Type-1 sum of squares and respective degrees of freedom.
#' 						In case of "sweep" funtion \code{\link{getSSQsweep}} will be called, otherwise, function \code{\link{getSSQqf}}                      
#' @param VarVC.method	(character) string specifying whether to use the algorithm given in Searle et al. (1992) which corresponds to \code{VarVC.method="scm"} or in
#' 						Giesbrecht and Burns (1985) which can be specified via "gb". Method "scm" (Searle, Casella, McCulloch)
#'                      is the exact algorithm but slower, "gb" (Giesbrecht, Burns) is termed "rough approximation"
#' 						by the authors, but sufficiently exact compared to e.g. SAS PROC MIXED (method=type1) which
#' 						uses the inverse of the Fisher-Information matrix as approximation. For balanced designs all
#'                      methods give identical results, only in unbalanced designs differences occur. 
#' @param MME			(logical) TRUE = (M)ixed (M)odel (E)quations will be solved, i.e. 'VCA' object will have additional elements
#' 						"RandomEffects", "FixedEffects", "VarFixed" (variance-covariance matrix of fixed effects) and the "Matrices"
#' 						element has addional elements corresponding to intermediate results of solving MMEs.
#' 						FALSE = do not solve MMEs, which reduces the computation time for very complex models significantly.
#' @param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' 
#' @return (object) of class 'VCA'
#' 
#' @aliases anovaVCA
#' 
#' @seealso \code{\link{anovaMM}}, \code{\link{remlVCA}}, \code{\link{remlMM}}, \code{\link{print.VCA}}, \code{\link{VCAinference}}, 
#' 			\code{\link{getCmatrix}}, \code{\link{ranef}}, \code{\link{plotRandVar}}, \code{\link{stepwiseVCA}}
#' 
#' @references 
#' 
#' Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York
#' 
#' Goodnight, J.H. (1979), A Tutorial on the SWEEP Operator, The American Statistician, 33:3, 149-158
#' 
#' Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#' Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486
#' 
#' Satterthwaite, F.E. (1946),  An Approximate Distribution of Estimates of Variance Components., 
#' Biometrics Bulletin 2, 110-114
#' 
#' Gaylor,D.W., Lucas,H.L., Anderson,R.L. (1970), Calculation of Expected Mean Squares by the Abbreviated Doolittle and Square Root Methods., 
#' Biometrics 26 (4): 641-655
#' 
#' SAS Help and Documentation PROC MIXED, SAS Institute Inc., Cary, NC, USA
#' 
#' SAS Help and Documentation PROC VARCOMP, SAS Institute Inc., Cary, NC, USA
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' \dontrun{
#' 
#' # load data (CLSI EP05-A2 Within-Lab Precision Experiment) 
#' data(dataEP05A2_2)
#' 
#' # perform ANOVA-estimation of variance components
#' res <- anovaVCA(y~day/run, dataEP05A2_2)
#' res
#' 
#' # desing with two main effects (ignoring the hierarchical structure of the design)
#' anovaVCA(y~day+run, dataEP05A2_2)
#' 
#' # compute confidence intervals, perform F- and Chi-Squared tests
#' INF <- VCAinference(res, total.claim=3.5, error.claim=2)
#' INF
#' 
#' ### load data from package
#' data(VCAdata1)
#' 
#' data_sample1 <- VCAdata1[VCAdata1$sample==1,]
#' 
#' ### plot data for visual inspection (there is no variance between runs on a day)
#' varPlot(y~lot/day/run, data_sample1)
#' 
#' ### estimate VCs for 4-level hierarchical design (error counted) for sample_1 data
#' anovaVCA(y~lot/day/run, data_sample1)
#' 
#' ### using different model (ignoring the hierarchical structure of the design)
#' anovaVCA(y~lot+day+lot:day:run, data_sample1)
#' 
#' ### same model with unbalanced data
#' anovaVCA(y~lot+day+lot:day:run, data_sample1[-c(1,11,15),])
#' 
#' ### use the numerical example from the CLSI EP05-A2 guideline (p.25)
#' data(Glucose)
#' res.ex <- anovaVCA(result~day/run, Glucose)
#' 
#' ### also perform Chi-Squared tests
#' ### Note: in guideline claimed SD-values are used, here, claimed variances are used
#' VCAinference(res.ex, total.claim=3.4^2, error.claim=2.5^2)
#' 
#' ### now use the six sample reproducibility data from CLSI EP5-A3
#' ### and fit per sample reproducibility model
#' data(CA19_9)
#' fit.all <- anovaVCA(result~site/day, CA19_9, by="sample")
#' 
#' reproMat <- data.frame(
#'  Sample=c("P1", "P2", "Q3", "Q4", "P5", "Q6"),
#'  Mean= c(fit.all[[1]]$Mean, fit.all[[2]]$Mean, fit.all[[3]]$Mean, 
#' 	        fit.all[[4]]$Mean, fit.all[[5]]$Mean, fit.all[[6]]$Mean),
#' 	Rep_SD=c(fit.all[[1]]$aov.tab["error","SD"], fit.all[[2]]$aov.tab["error","SD"],
#' 	         fit.all[[3]]$aov.tab["error","SD"], fit.all[[4]]$aov.tab["error","SD"],
#'           fit.all[[5]]$aov.tab["error","SD"], fit.all[[6]]$aov.tab["error","SD"]),
#' 	Rep_CV=c(fit.all[[1]]$aov.tab["error","CV[%]"],fit.all[[2]]$aov.tab["error","CV[%]"],
#'           fit.all[[3]]$aov.tab["error","CV[%]"],fit.all[[4]]$aov.tab["error","CV[%]"],
#'           fit.all[[5]]$aov.tab["error","CV[%]"],fit.all[[6]]$aov.tab["error","CV[%]"]),
#'  WLP_SD=c(sqrt(sum(fit.all[[1]]$aov.tab[3:4,"VC"])),sqrt(sum(fit.all[[2]]$aov.tab[3:4, "VC"])),
#'           sqrt(sum(fit.all[[3]]$aov.tab[3:4,"VC"])),sqrt(sum(fit.all[[4]]$aov.tab[3:4, "VC"])),
#'           sqrt(sum(fit.all[[5]]$aov.tab[3:4,"VC"])),sqrt(sum(fit.all[[6]]$aov.tab[3:4, "VC"]))),
#'  WLP_CV=c(sqrt(sum(fit.all[[1]]$aov.tab[3:4,"VC"]))/fit.all[[1]]$Mean*100,
#'           sqrt(sum(fit.all[[2]]$aov.tab[3:4,"VC"]))/fit.all[[2]]$Mean*100,
#'           sqrt(sum(fit.all[[3]]$aov.tab[3:4,"VC"]))/fit.all[[3]]$Mean*100,
#'           sqrt(sum(fit.all[[4]]$aov.tab[3:4,"VC"]))/fit.all[[4]]$Mean*100,
#'           sqrt(sum(fit.all[[5]]$aov.tab[3:4,"VC"]))/fit.all[[5]]$Mean*100,
#'           sqrt(sum(fit.all[[6]]$aov.tab[3:4,"VC"]))/fit.all[[6]]$Mean*100),
#'  Repro_SD=c(fit.all[[1]]$aov.tab["total","SD"],fit.all[[2]]$aov.tab["total","SD"],
#'             fit.all[[3]]$aov.tab["total","SD"],fit.all[[4]]$aov.tab["total","SD"],
#'             fit.all[[5]]$aov.tab["total","SD"],fit.all[[6]]$aov.tab["total","SD"]),
#'  Repro_CV=c(fit.all[[1]]$aov.tab["total","CV[%]"],fit.all[[2]]$aov.tab["total","CV[%]"],
#'             fit.all[[3]]$aov.tab["total","CV[%]"],fit.all[[4]]$aov.tab["total","CV[%]"],
#'             fit.all[[5]]$aov.tab["total","CV[%]"],fit.all[[6]]$aov.tab["total","CV[%]"]))
#'  
#'  for(i in 3:8) reproMat[,i] <- round(reproMat[,i],digits=ifelse(i%%2==0,1,3))
#'  reproMat
#' 
#' # now plot the precision profile over all samples
#' plot(reproMat[,"Mean"], reproMat[,"Rep_CV"], type="l", main="Precision Profile CA19-9",
#' 		xlab="Mean CA19-9 Value", ylab="CV[%]")
#' grid()
#' points(reproMat[,"Mean"], reproMat[,"Rep_CV"], pch=16)
#' 
#' # load another example dataset and extract the "sample==1" subset
#' data(VCAdata1)
#' sample1 <- VCAdata1[which(VCAdata1$sample==1),]
#'  
#' # generate an additional factor variable and random errors according to its levels
#' sample1$device <- gl(3,28,252)                                      
#' set.seed(505)
#' sample1$y <- sample1$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     
#' 
#' # fit a crossed-nested design with main factors 'lot' and 'device' 
#' # and nested factors 'day' and 'run' nested below 
#' res1 <- anovaVCA(y~(lot+device)/day/run, sample1) 
#' res1
#' 
#' # fit same model for each sample using by-processing
#' lst <- anovaVCA(y~(lot+device)/day/run, VCAdata1, by="sample")
#' lst
#' 
#' # now fitting a nonsense model on the complete dataset "VCAdata1" 
#' # using both methods for fitting ANOVA Type-1 sum of squares
#' # SSQ.method="qf" used to be the default up to package version 1.1.1
#' # took ~165s on a Intel Xeon E5-2687W (3.1GHz) in V1.1.1, now takes ~95s
#' system.time(res.qf <- anovaVCA(y~(sample+lot+device)/day/run, VCAdata1, SSQ.method="qf"))
#' # the SWEEP-operator is the new default since package version 1.2
#' # takes ~5s
#' system.time(res.sw <- anovaVCA(y~(sample+lot+device)/day/run, VCAdata1, SSQ.method="sweep"))
#' # applying functions 'anova' and 'lm' in the same manner takes ~ 265s
#' system.time(res.lm <- anova(lm(y~(sample+lot+device)/day/run, VCAdata1)))
#' res.qf
#' res.sw
#' res.lm
#' }

anovaVCA <- function(	form, Data, by=NULL, NegVC=FALSE, SSQ.method=c("sweep", "qf"), 
		VarVC.method=c("gb", "scm"), MME=FALSE, quiet=FALSE)
{
	if(!is.null(by))
	{
		stopifnot(is.character(by))
		stopifnot(by %in% colnames(Data))
		stopifnot(is.factor(by) || is.character(by))
		
		levels  <- unique(Data[,by])
		res <- lapply(levels, function(x) anovaVCA(form=form, Data[Data[,by] == x,], NegVC=NegVC, SSQ.method=SSQ.method, VarVC.method=VarVC.method, MME=MME, quiet=quiet))
		names(res) <- paste(by, levels, sep=".")
		return(res)
	}
	
	stopifnot(class(form) == "formula")
	stopifnot(is.data.frame(Data))
	stopifnot(nrow(Data) > 2)                                               # at least 2 observations for estimating a variance
	stopifnot(is.logical(NegVC))
	
	VarVC.method <- match.arg(VarVC.method)
	SSQ.method   <- match.arg(SSQ.method)
	VarVC.method <- ifelse(SSQ.method == "sweep", "gb", VarVC.method)		# always use "gb", since A-matrices will not be computed	
	tobj <- terms(form, simplify=TRUE, keep.order=TRUE)                     # expand nested factors if necessary, retain ordering of terms in the formula
	form <- formula(tobj)
	
	res       	<- list()
	res$call  	<- match.call()
	res$Type  	<- "Random Model"
	res$data  	<- Data
	res$terms 	<- tobj
	res$VCnames <- c(attr(tobj, "term.labels"), "error")
	
	int <- res$intercept <- attr(tobj, "intercept") == 1	
	
	if(!attr(tobj, "response"))
		stop("You need to include a response variable in the formula!")
	resp <- as.character(form)[2]
	res$response <- resp
	
	stopifnot(resp %in% colnames(Data))
	stopifnot(is.numeric(Data[,resp]))
	
	Ndata <- nrow(Data)
	
	rmInd <- integer()
	
	resp.NA <- is.na(Data[,resp])
	
	if(any(resp.NA))
	{    
		rmInd <- c(rmInd, which(resp.NA))
		if(!quiet)
			warning("There are ", length(which(resp.NA))," missing values for the response variable (obs: ", paste(which(resp.NA), collapse=", "), ")!")
	}    
	
	fac  <- attr(tobj, "term.labels")
	vars <- rownames(attr(tobj, "factors"))[-1]                             # remove response
	Nvc  <- length(fac) + 1 
	
	if(!is.null(vars))
	{
		for(i in vars)                                                          # convert all nested factors as factor-objects
		{
			if( any(is.na(Data[,i])))
			{
				NAind <- which(is.na(Data[,i]))
				rmInd <- c(rmInd, NAind)
				if(!quiet)
					warning("Variable '", i,"' has ",length(NAind)," missing values (obs: ", paste(NAind, collapse=", "), ")!" )
			}
			if(length(levels(Data[,i])) >= 0.75*nrow(Data) && !quiet)
				warning("Variable >>> ", i," <<< has at least 0.75 * nrow(Data) levels!")
		}
		rmInd <- unique(rmInd)
	}	
	
	if(length(rmInd) > 0)
		Data <- Data[-rmInd,]
	
	for(i in rev(vars))														# order data 														
	{
		tmp <- Data[,i]
		tmp.num <- suppressWarnings(as.numeric(as.character(tmp)))			# warnings are very likely here
		
		if(!any(is.na(tmp.num)))											# all levels are numbers
		{
			Width 	 <- max(nchar(as.character(unique(Data[,i])))) + 1
			tmp   	 <- suppressWarnings(formatC(Data[,i], width=Width))
			uLvl  	 <- unique(Data[,i])
			Data  	 <- Data[order(tmp),]
			Data[,i] <- factor(as.character(Data[,i]), levels=uLvl[order(unique(tmp))])
		}		
		else
			Data[,i] <- factor(Data[,i])									# automatically orders
	}

	Mean <- mean(Data[,resp], na.rm=TRUE)                                   # mean computed after removing incomplete observations
	vcol <- rownames(attr(tobj, "factors"))
	if(is.null(vcol))
		vcol <- resp
	Data <- na.omit(Data[,vcol, drop=F])
	Nobs <- N <- nrow(Data)
	
	if(SSQ.method == "qf")
		tmp.res <- getSSQqf(Data, tobj)										# determine ANOVA Type-1 sum of squares using the quadratic form approach									
	else
		tmp.res <- getSSQsweep(Data, tobj)									# determine ANOVA Type-1 sum of squares using sweeping
	
	gc(verbose=FALSE)
	
	Lmat    <- tmp.res$Lmat
	aov.tab <- tmp.res$aov.tab												# basic ANOVA-table
	DF <- aov.tab[,"DF"]
	SS <- aov.tab[,"SS"]
	C <- getCmatrix(form, Data, aov.tab[,"DF"], "SS", MM=Lmat$Zre) 			# compute coefficient matrix C in ss = C * s
	Ci  <- solve(C)
	C2  <- apply(C, 2, function(x) x/DF)                                    # coefficient matrix for mean squares (MS)
	Ci2 <- solve(C2)
	VC  <- as.matrix( Ci %*% SS)                                            # solve for VCs (p.173)
	aov.tab$VC <- VC
	colnames(aov.tab) <- c("DF", "SS", "MS", "VC")
	
	IndNegVC <- 1:nrow(aov.tab)                                             # should negative VC-estimates contribute to total variance 
	NegVCmsg <- ""                                                          # message text in case that VC is negative
	VCoriginal <- aov.tab[,"VC"]                                            # keep original ANOVA-estimates of VCs
	if(!NegVC)
	{
		IndNegVC <- which(aov.tab[,"VC"] < 0)  
		if(length(IndNegVC) > 0)                                            # there are negative VC
		{
			aov.tab[IndNegVC, "VC"] <- 0                                    # set negative VC to 0
			NegVCmsg <- "* VC set to 0"
		}
	}    
	totVC  <- sum(aov.tab$VC)
	
	aov.tab <- rbind(total=c(NA, NA, NA, totVC), aov.tab)    
	
	aov.tab["total", "DF"] <- SattDF(c(C2 %*% aov.tab[-1,"VC"]), Ci=Ci2, DF=DF)   	# will automatically adapt ANOVA-MS if any VCs were set to 0 
	
	suppressWarnings(aov.tab <- cbind(aov.tab, SD=sqrt(aov.tab[,"VC"])))    			# warnings suppressed because sqrt of negative numbers doese not exists
	aov.tab <- cbind(aov.tab, "CV[%]"=aov.tab[,"SD"]*100/Mean)
	aov.tab <- cbind(aov.tab, "%Total"=aov.tab[,"VC"]*100/totVC)
	aov.tab <- aov.tab[,c("DF", "SS", "MS", "VC", "%Total", "SD", "CV[%]")]    
	aov.tab <- as.matrix(aov.tab)
	aov.tab <- apply(aov.tab, 1:2, function(x) ifelse(is.nan(x), NA, x))
	
	res$aov.tab <- aov.tab
	res$rmInd   <- rmInd
	
	res$Nvc			 <- Nvc
	res$VarVC.method <- VarVC.method
	res$SSQ.method   <- SSQ.method
	res$Mean         <- Mean
	res$Nobs         <- Nobs
	res$EstMethod    <- "ANOVA"
	res$VCoriginal   <- VCoriginal
	res$NegVCmsg     <- NegVCmsg
	res$NegVC        <- NegVC
	if(Nobs != Ndata)
		res$Nrm <- Ndata - Nobs                                	# save number of observations that were removed due to missing data
	
	Lmat$C.SS	<- C
	Lmat$C.MS	<- C2
	Lmat$Ci.SS  <- Ci
	Lmat$Ci.MS  <- Ci2
	Lmat$y      <- matrix(Data[, resp], ncol=1)								# column vector of observations
	Lmat$X      <- matrix(ifelse(int, 1, 0), ncol=1, nrow=nrow(Data))		# design matrix of fixed effects
	colnames(Lmat$X) <- "int"
	if(int)
	{
		res$fe.assign <- 0
		attr(res$fe.assign, "terms") <- "int"
	}
	Lmat$b	    <- matrix(Mean, 1, 1)										# column vector of fixed effects (intercept here for random models)
	Lmat$rf.ind <- 1:(nrow(aov.tab)-1)
	Lmat$VCall  <- aov.tab[-1, "VC"]
	res$Matrices <- Lmat                                   	# needed for variance-covariance matrix of VCs
	
	res$formula  <- form
	res$balanced <- if(isBalanced(form, Data)) 
				"balanced"  
			else 
				"unbalanced"
	class(res) <- "VCA"    
	
	if(MME)
		res <- solveMME(res)
	
	gc(verbose=FALSE)													# trigger garbage collection
	return(res)
}



#' Extract Residuals of a 'VCA' Object.
#' 
#' Function extracts marginal or conditional residuals from a 'VCA' object, 
#' representing a linear mixed model.
#' 
#' There are two types of residuals which can be extraced from a 'VCA' object.
#' Marginal residuals correspond to \eqn{e_m = y - \hat{y}}{e_m = y - y_hat}, where \eqn{\hat{y} = Xb}{y_hat = Xb} with \eqn{X}
#' being the design matrix of fixed effects and \eqn{b} being the column vector of fixed
#' effects parameter estimates. Conditional residuals are defined as \eqn{e_c = y - Xb - Zg},
#' where \eqn{Z} corresponds to the designs matrix of random effects \eqn{g}. 
#' Whenever 'obj' is a pure-error model, e.g. 'y~1' both options will return the same values 
#' \eqn{y - Xb} and \eqn{b} corresponds to the intercept.
#' Each type of residuals can be standardized, studentized, or transformed to pearson-type residuals. 
#' The former corresponds to a transformation of residuals to have mean 0 and variance equal to 1 (\eqn{(r - \bar{r})/\sigma_{r}}{[r - mean(r)]/sd(r)]}). 
#' Studentized residuals emerge from dividing raw residuals by the square-root of diagonal elements of the corresponding 
#' variance-covariance matrix. For conditional residuals, this is \eqn{Var(c) = P = RQR}, with \eqn{Q = V^{-1}(I - H)}{Q = V"(I - H)},
#' \eqn{H = XT} being the hat-matrix, and \eqn{T = (X^{T}V^{-1}X)^{-1}X^{T}V^{-1}}{T = (X'V"X)"X'V"}. For marginal residuals, this matrix
#' is \eqn{Var(m) = O = V - Q}. Here, >\eqn{^{T}}{'}< denotes the matrix transpose operator, 
#' and >\eqn{^{-1}}{"}< the regular matrix inverse. Pearson-type residuals are computed in the same manner as studentized, only
#' the variance-covariance matrices differ. For marginal residuals this is equal to \eqn{Var(y) = V}, for conditional residuals
#' this is \eqn{Var(c) = R} (see \code{\link{getV}} for details). 
#' 
#' @param object		(VCA) object
#' @param type			(character) string specifying the type of residuals to be returned,
#'                      valid options are "marginal" and "conditional" or abbreviations
#' @param mode			(character) string or abbreviation specifying the specific transformation
#'                      applied to a certain type of residuals. There are "raw" (untransformed), 
#'                      "standardized", "studentized" and "pearson" (see details) residuals.
#' @param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#' @param ...			additional parameters
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @method residuals VCA
#' @S3method residuals VCA
#' 
#' @references 
#' 
#' Hilden-Minton, J. A. (1995). Multilevel diagnostics for mixed and hierarchical linear
#' models. Dissertation, University of California, Los Angeles.
#' 
#' Nobre, J. S. & Singer, J. M. (2007). Residual analysis for linear mixed models. Biometrical
#' Journal, 49, 863-875.
#' 
#' Schuetzenmeister, A. and Piepho, H.P. (2012). Residual analysis of linear mixed models using a simulation approach.
#' Computational Statistics and Data Analysis, 56, 1405-1416
#' 
#' @examples
#' \dontrun{
#' data(VCAdata1)
#' datS1 <- VCAdata1[VCAdata1$sample==1,]
#' fit1  <- anovaVCA(y~(lot+device)/(day)/(run), datS1) 
#' 
#' # default is conditional (raw) residuals
#' resid(fit1)
#' resid(fit1, "m")
#' 
#' # get standardized version
#' resid(fit1, mode="stand")		# conditional residuals (default)
#' resid(fit1, "marg", "stand")		# marginal residuals
#' 
#' # get studentized version, taking their 
#' # covariances into account
#' resid(fit1, mode="stud")		# conditional residuals (default)
#' resid(fit1, "marg", "stud")		# marginal residuals
#' }
#' 
#' @aliases resid
#' 
#' @seealso \code{\link{ranef}}, \code{\link{anovaVCA}}, \code{\link{anovaMM}}

residuals.VCA <- function(object, type=c("conditional", "marginal"), mode=c("raw", "student", "standard", "pearson"), quiet=FALSE, ...)
{		
	Call <- match.call()
	obj <- object
	
	if(is.list(obj) && class(obj) != "VCA")
	{
		if(!all(sapply(obj, class) == "VCA"))
			stop("Only lists of 'VCA' object are accepted!")
		
		obj.len <- length(obj)
		
		if(!"msgEnv" %in% ls(.GlobalEnv))
			msgEnv <<- new.env(parent=emptyenv())
		
		assign("VCAinference.obj.is.list", TRUE, envir=msgEnv)			# indicate that a list-type object was passed intially
				
		res <- mapply(	FUN=residuals.VCA, object=obj,
				type=type[1], mode=mode[1], SIMPLIFY=FALSE)
		
		names(res) <- names(obj)
		
		if(obj.len == 1)			# mapply returns a list of length 2 in case that length(obj) was equal to 1
			res <- res[1]
		
		rm("VCAinference.obj.is.list", envir=msgEnv)
		
		return(res)
	}	
	
	stopifnot(class(obj) == "VCA")
	
	type <- match.arg(type)
	mode <- match.arg(mode)
	
	if(is.null(obj$FixedEffects))
	{
		obj  <- solveMME(obj)
		nam0 <- deparse(Call$object)
		nam1 <- sub("\\[.*", "", nam0)
		
		if(length(nam1) == 1 && nam1 %in% names(as.list(.GlobalEnv)))
		{
			expr <- paste(nam0, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else
		{
			if( !"VCAinference.obj.is.list" %in% names(as.list(msgEnv)) && !quiet)
				warning("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
		}
	}
	
	Xb <- getMat(obj, "X") %*% obj$FixedEffects
	y  <- getMat(obj, "y")
	
	if(as.character(obj$terms)[3] == "1")							# pure error model
	{
		type <- "marginal"
	}
	
	if(type == "marginal")											# marginal residuals
	{
		res <- y - Xb
		V   <- getMat(obj, "V")
		
		if(mode == "student")
		{			
			X  <- getMat(obj, "X")
			Vi <- getMat(obj, "Vi")
			
			#Q  <- X %*% MPinv(t(X) %*% Vi %*% X) %*% t(X)
			K   <- obj$VarFixed																	
			Q   <- X %*% K %*% t(X)
			res <- res/sqrt(diag(V-Q))								# apply studentization
		}
		else if(mode == "pearson")
		{
			res <- res/sqrt(diag(V))								# Pearson-type residuals
		}
	}
	else															# conditional residuals
	{
		Z   <- getMat(obj, "Z")
		if(obj$EstMethod == "REML")									# re-order Z to align with random effects
			Z <- Z[,unlist(obj$ColOrderZ)]
		
		Zg  <- Z %*% ranef(obj)										# order according to column-order of Z
		
		res <- y - Xb - Zg
		
		R  	<- getMat(obj, "R")
		
		if(mode == "student")
		{
			mats <- obj$Matrices
			
			Q	 <- mats$Q
			
			if(is.null(Q))
			{
				X  <- mats$X
				T  <- mats$T
				Vi <- mats$Vi
				mats$H <- H  <- X %*% T
				mats$Q <- Q  <- Vi %*% (diag(nrow(H))-H)
				
				nam0 <- deparse(Call$object)	
				nam1 <- sub("\\[.*", "", nam0) 
				
				if(length(nam1) == 1 && nam1 %in% names(as.list(.GlobalEnv)))	# write back to object in the calling env
				{
					obj$Matrices <- mats
					expr <- paste(nam0, "<<- obj")
					eval(parse(text=expr))
				}
			}
			P 	 <- R %*% Q %*% R
			res	 <- res / sqrt(diag(P))								# apply studentization
		}
		else if(mode == "pearson")
		{
			res <- res/sqrt(diag(R))								# Pearson-type residuals
		}
	}
	res <- c(res[,1])											
	
	if(mode=="standard")
	{
		res <- res / sd(res)										# apply standardization
	}
	names(res) <- rownames(obj$data)		
	attr(res, "type") <- type
	attr(res, "mode") <- mode
	return(res)
}






#' Bottom-Up Step-Wise VCA-Analysis of the Complete Dataset.
#' 
#' Function performs step-wise VCA-analysis on a fitted VCA-object by leaving out N-1 to 0
#' top-level variance components (VC).
#' 
#' This function uses the complete data to quantify sub-sets of variance components.
#' In each step the current total variance is estimated by subtracting the sum of all left-out VCs
#' from the total variance of the initial VCA object. Doing this guarantees that the contribution to the total 
#' variance which is due to left-out VCs is accounted for, i.e. it is estimated but not included/reported.
#' The degrees of freedom (DFs) the emerging total variances of sub-sets are estimated using the Satterthwaite
#' approximation. This is achieved by extracting the corresponding sub-matrix from the coefficient matrix \eqn{C} of
#' the 'VCA' object, the sub-vector of ANOVA mean squares, and the sub-vector of degrees of freedom and calling
#' function \code{\link{SattDF}} method="total".
#' This step-wise procedure starts one-level above error (repeatability) and ends at the level of the upper-most VC.
#' 
#' @param obj			(VCA) object representing the complete analysis
#' @param VarVC			(logical) TRUE = estimate complete variance-covariance matrix of variance components and include
#'                      corresponding sub-matrices in result-objects of the step-wise analysis. This allows to cosntruct
#'                      confidence intervals on all VCs via function \code{\link{VCAinference}}
#' @param VarVC.method	(character) string specifying the algorithm to be used for estimating variance-covariance matrix
#'                      of VCs (see \code{\link{anovaMM}} for details).
#' 
#' @return (list) of (simplified) 'VCA' objects representing analysis-result of sub-models
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' \dontrun{
#' data(VCAdata1)
#' datS7L1 <- VCAdata1[VCAdata1$sample == 7 & VCAdata1$lot == 1, ]
#' fit0 <- anovaVCA(y~device/day/run, datS7L1, MME=TRUE)
#' 
#' # complete VCA-analysis result
#' fit0
#' 
#' # perform step-wise (bottom-up) VCA-analyses
#' sw.res <- stepwiseVCA(fit0, VarVC=TRUE)
#' sw.res
#' 
#' # get CIs on intermediate precision 
#' VCAinference(sw.res[["device:day"]], VarVC=TRUE)
#' }

stepwiseVCA <- function(obj, VarVC=FALSE, VarVC.method=c("scm", "gb"))
{
	stopifnot(obj$Type == "Random Model")					# only random models correspond to true VCA-analyses
	VarVC.method <- match.arg(VarVC.method)
	VarVC.method <- ifelse(obj$SSQ.method == "sweep", "gb", VarVC.method)		# always use "gb", since A-matrices will not be computed	
	tab <- obj$aov.tab
	Nstp  <- nrow(tab)-2
	res <- vector("list", length=Nstp)
	names(res) <- rev(rownames(tab)[2:(Nstp+1)])
	Ci <- getMat(obj, "Ci.MS")
	N  <- nrow(tab)
	if(VarVC)
		VarCov <- vcovVC(obj, method=VarVC.method)
	for(i in 1:Nstp)
	{
		C <- Ci[(N-i-1):(N-1), (N-i-1):(N-1)]
		M <- tab[(N-i):(N), "MS"]
		D <- tab[(N-i):(N), "DF"]
		
		DF <- SattDF(M, C, D)
		mat <- tab[(N-i):N, ]	
		mat <- rbind(total=c(DF=DF, SS=NA, MS=NA, VC=sum(mat[,"VC"]), "%Total"=100, SD=NA, "CV[%]"=NA), mat)
		mat[, "%Total"] <- 100 * mat[,"VC"] / mat["total", "VC"]
		mat["total", "SD"] <- sqrt(mat["total", "VC"])
		mat["total", "CV[%]"] <- 100 * mat["total", "SD"] / obj$Mean
		ele <- list(aov.tab=mat)
		ele$NegVCmsg <- obj$NegVCmsg
		class(ele) <- "VCA"
		ele$Type <- "Random Model"
		ele$Mean <- obj$Mean
		ele$Nobs <- obj$Nobs
		ele$balanced <- obj$balanced
		ele$Matrices <- list(Ci.MS=C, rf.ind=1:(nrow(mat)-1))
		if(VarVC)											# extract sub-matrix of the variance-covariance matrix of VCs
		{
			ele$VarCov <- VarCov[(N-i-1):(N-1), (N-i-1):(N-1)]
			ele$VarVC.method <- VarVC.method
		}
		
		res[[i]] <- ele
	}
	return(res)
}
