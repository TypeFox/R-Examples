### rirls.spls.R  (2014-10)
###
###    Ridge Iteratively Reweighted Least Squares followed by Adaptive Sparse PLS regression for binary response
###
### Copyright 2014-10 Ghislain DURIF
###
### Adapted from rpls function in plsgenomics package, copyright 2006-01 Sophie Lambert-Lacroix
###
### This file is part of the `plsgenomics' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


rirls.spls <- function(Xtrain, Ytrain, lambda.ridge, lambda.l1, ncomp, Xtest=NULL, adapt=TRUE, maxIter=100, svd.decompose=TRUE) {
	
	
	#####################################################################
	#### Initialisation
	#####################################################################
	Xtrain <- as.matrix(Xtrain)
	ntrain <- nrow(Xtrain) # nb observations
	p <- ncol(Xtrain) # nb covariates
	index.p <- c(1:p)
	Ytrain <- as.matrix(Ytrain)
	q = ncol(Ytrain)
	one <- matrix(1,nrow=1,ncol=ntrain)
	
	
	#####################################################################
	#### Tests on type input
	#####################################################################
	# On Xtrain
	if ((!is.matrix(Xtrain)) || (!is.numeric(Xtrain))) {
		stop("Message from rirls.spls: Xtrain is not of valid type")
	}
	
	if (p==1) {
		stop("Message from rirls.spls: p=1 is not valid")}
	
	# On Xtest if necessary
	if (!is.null(Xtest)) {
		
		if (is.vector(Xtest)==TRUE) {
			Xtest <- matrix(Xtest,nrow=1)
		}
		
		Xtest <- as.matrix(Xtest)
		ntest <- nrow(Xtest) 
		
		if ((!is.matrix(Xtest)) || (!is.numeric(Xtest))) {
			stop("Message from rirls.spls: Xtest is not of valid type")}
		
		if (p != ncol(Xtest)) {
			stop("Message from rirls.spls: columns of Xtest and columns of Xtrain must be equal")
		}	
	}
	
	# On Ytrain
	if ((!is.matrix(Ytrain)) || (!is.numeric(Ytrain))) {
		stop("Message from rirls.spls: Ytrain is not of valid type")
	}
	
	if (q != 1) {
		stop("Message from rirls.spls: Ytrain must be univariate")
	}
	
	if (nrow(Ytrain)!=ntrain) {
		stop("Message from rirls.spls: the number of observations in Ytrain is not equal to the Xtrain row number")
	}
	
	# On Ytrain value
	if (sum(is.na(Ytrain))!=0) {
		stop("Message from rirls.spls: NA values in Ytrain")
	}
	
	if (sum(!(Ytrain %in% c(0,1)))!=0) {
		stop("Message from rirls.spls: Ytrain is not of valid type")
	}
	
	if (sum(as.numeric(table(Ytrain))==0)!=0) {
		stop("Message from rirls.spls: there are empty classes")
	}
	
	# On hyper parameter: lambda.ridge, lambda.l1
	if ((!is.numeric(lambda.ridge)) || (lambda.ridge<0) || (!is.numeric(lambda.l1)) || (lambda.l1<0)) {
		stop("Message from rirls.spls: lambda is not of valid type")
	}
	
	# ncomp type
	if ((!is.numeric(ncomp)) || (round(ncomp)-ncomp!=0) || (ncomp<0) || (ncomp>p)) {
		stop("Message from rirls.spls: ncomp is not of valid type")
	}
	
	# maxIter
	if ((!is.numeric(maxIter)) || (round(maxIter)-maxIter!=0) || (maxIter<1)) {
		stop("Message from rirls.spls: maxIter is not of valid type")
	}
	
	
	#####################################################################
	#### Move into the reduced space
	#####################################################################
	
	r <- min(p, ntrain)
	DeletedCol <- NULL
	
	### Standardize the Xtrain matrix
	# standard deviation (biased one) of Xtrain
	sigma2train <- apply(Xtrain, 2, var) * (ntrain-1)/(ntrain)
	
	# test on sigma2train
	# predictor with null variance ?
	if (sum(sigma2train < .Machine$double.eps)!=0){
		
		# predicteur with non null variance < 2 ?
		if (sum(sigma2train < .Machine$double.eps)>(p-2)){
			stop("Message from rirls.spls: the procedure stops because number of predictor variables with no null variance is less than 1.")
		}
		
		warning("There are covariables with nul variance")
		
		# remove predictor with null variance
		Xtrain <- Xtrain[,which(sigma2train >= .Machine$double.eps)]
		if (!is.null(Xtest)) {
			Xtest <- Xtest[,which(sigma2train>= .Machine$double.eps)]
		}
		
		# list of removed predictors
		DeletedCol <- index.p[which(sigma2train < .Machine$double.eps)]
		
		# removed null standard deviation
		sigma2train <- sigma2train[which(sigma2train >= .Machine$double.eps)]
		
		# new number of predictor
		p <- ncol(Xtrain)
		r <- min(p,ntrain)
	}
	
	# mean of Xtrain
	meanXtrain <- apply(Xtrain,2,mean)
	
	# center and scale Xtrain
	sXtrain <- scale(Xtrain, center=meanXtrain, scale=sqrt(sigma2train))
	
	sXtrain.nosvd = sXtrain # keep in memory if svd decomposition
	
	# Compute the svd when necessary -> case p > ntrain (high dim)
	if ((p > ntrain) && (svd.decompose)) {
		# svd de sXtrain
		svd.sXtrain <- svd(t(sXtrain))
		# number of singular value non null
		r <- length(svd.sXtrain$d[abs(svd.sXtrain$d)>10^(-13)])
		V <- svd.sXtrain$u[,1:r]
		D <- diag(c(svd.sXtrain$d[1:r]))
		U <- svd.sXtrain$v[,1:r]
		sXtrain <- U %*% D
		rm(D)
		rm(U)
		rm(svd.sXtrain)
	}
	
	# center and scale Xtest if necessary
	if (!is.null(Xtest)) {
		
		meanXtest <- apply(Xtest,2,mean)
		sigma2test <- apply(Xtest,2,var)
		
		sXtest <- scale(Xtest, center=meanXtrain, scale=sqrt(sigma2train))
		
		sXtest.nosvd <- sXtest # keep in memory if svd decomposition
		
		# if svd decomposition
		if ((p > ntrain) && (svd.decompose)) {
			sXtest <- sXtest%*%V
		}
	}
	
	
	#####################################################################
	#### Ridge IRLS step
	#####################################################################
	
	fit <- wirrls(Y=Ytrain, Z=cbind(rep(1,ntrain),sXtrain), Lambda=lambda.ridge, NbrIterMax=maxIter, WKernel=diag(rep(1,ntrain)))
	
	converged=fit$Cvg
	
	#  Check WIRRLS convergence
	if (converged==0) {
		warning("Message from rirls.spls : Ridge IRLS did not converge; try another lambda.ridge value")
	}
	
	# if ncomp == 0 then wirrls without spls step
	if (ncomp==0) {
		BETA <- fit$Coefficients
	}
	
	
	#####################################################################
	#### weighted SPLS step
	#####################################################################
	
	# if ncomp > 0
	if (ncomp!=0) {
		
		#Compute ponderation matrix V and pseudo variable z
		#Pseudovar = Eta + W^-1 Psi
		
		# Eta = X * betahat (covariate summary)
		Eta <- cbind(rep(1, ntrain), sXtrain) %*% fit$Coefficients
		
		## Run SPLS on Xtrain without svd decomposition
		sXtrain = sXtrain.nosvd
		if (!is.null(Xtest)) {
			sXtest = sXtest.nosvd
		}
		
		# mu = h(Eta)
		mu = 1 / (1 + exp(-Eta))
		
		# ponderation matrix : V
		diagV <- mu * (1-mu)
		V <- diag(c(diagV))
		
		# inv de V
		diagVinv = 1/ifelse(diagV!=0, diagV, diagV+0.00000001)
		Vinv = diag(c(diagVinv))
		
		Psi <- Ytrain-mu
		
		pseudoVar = Eta + Vinv %*% Psi
		
		# V-Center the sXtrain and pseudo variable
		sumV=sum(diagV)
		
		# Weighted centering of Pseudo variable
		VmeanPseudoVar <- sum(V %*% Eta + Psi ) / sumV
		
		# Weighted centering of sXtrain
		VmeansXtrain <- t(diagV)%*%sXtrain/sumV
		
		# SPLS(X, pseudo-var, weighting = V)
		resSPLS = spls.adapt(Xtrain=sXtrain, Ytrain=pseudoVar, ncomp=ncomp, weight.mat=V, lambda.l1=lambda.l1, adapt=adapt, center.X=TRUE, scale.X=FALSE, center.Y=TRUE, scale.Y=FALSE, weighted.center=TRUE)
		
		BETA = resSPLS$betahat.nc
		
	}
	
	
	#####################################################################
	#### classification step
	#####################################################################
	
	hatY <- numeric(ntrain)
	
	if (!is.null(Xtest)) {
		
		hatYtest <- cbind(rep(1,ntest),sXtest) %*% BETA
		hatYtest <- as.numeric(hatYtest>0)
		
		proba.test = inv.logit( cbind(rep(1,ntest),sXtest) %*% BETA )
		
		hatY <- cbind(rep(1,ntrain),sXtrain) %*% BETA
		hatY <- as.numeric(hatY>0)
		
	} else {
		
		hatYtest <- NULL
		
		proba.test <- NULL
		
		hatY <- cbind(rep(1,ntrain),sXtrain) %*% BETA
		hatY <- as.numeric(hatY>0)
		
	}
	
	
	#####################################################################
	#### Conclude
	#####################################################################
	
	Coefficients=BETA
	
	Coefficients[-1] <- diag(c(1/sqrt(sigma2train)))%*%BETA[-1]
	
	Coefficients[1] <- BETA[1] - meanXtrain %*% Coefficients[-1]
	
	
	#### RETURN
	
	result <- list(Coefficients=Coefficients, hatY=hatY, hatYtest=hatYtest, DeletedCol=DeletedCol, A=resSPLS$A, converged=converged, X.score=resSPLS$X.score, X.weight=resSPLS$X.weight, sXtrain=resSPLS$sXtrain, sPseudoVar=resSPLS$sYtrain, lambda.ridge=lambda.ridge, lambda.l1=lambda.l1, ncomp=ncomp, V=resSPLS$V, proba.test=proba.test, Xtrain=Xtrain, Ytrain=Ytrain)
	class(result) <- "rirls.spls"
	return(result)
	
	
	
	
}