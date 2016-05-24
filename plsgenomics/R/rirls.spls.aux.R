### rirls.spls.aux.R  (2014-10)
###
###    Ridge Iteratively Reweighted Least Squares followed by Adaptive Sparse PLS regression for binary responser
###    Short version for multiple call in cross-validation procedure
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


rirls.spls.aux <- function(sXtrain, sXtrain.nosvd=NULL, Ytrain, lambda.ridge, lambda.l1, ncomp, sXtest, sXtest.nosvd=NULL, adapt=TRUE, maxIter=100, svd.decompose=TRUE, meanXtrain, sigma2train) {
	
	
	#####################################################################
	#### Initialisation
	#####################################################################
	sXtrain <- as.matrix(sXtrain)
	ntrain <- nrow(sXtrain) # nb observations
	p <- ncol(sXtrain) # nb covariates
	index.p <- c(1:p)
	Ytrain <- as.matrix(Ytrain)
	q = ncol(Ytrain)
	one <- matrix(1,nrow=1,ncol=ntrain)
	ntest <- nrow(sXtest)
	
	
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
		sXtest = sXtest.nosvd
		
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
	
	hatYtest <- cbind(rep(1,ntest),sXtest) %*% BETA
	hatYtest <- as.numeric(hatYtest>0)
	
	
	#####################################################################
	#### Conclude
	#####################################################################
	
	Coefficients=BETA
	
	Coefficients[-1] <- diag(c(1/sqrt(sigma2train)))%*%BETA[-1]
	
	Coefficients[1] <- BETA[1] - meanXtrain %*% Coefficients[-1]
	
	
	#### RETURN
	
	result <- list(Coefficients=Coefficients, hatYtest=hatYtest, converged=converged)
	class(result) <- "rirls.spls.aux"
	return(result)
		
}
