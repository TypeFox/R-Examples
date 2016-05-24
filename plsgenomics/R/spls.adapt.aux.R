### spls.adapt.aux.R  (2014-10)
###
###    Adaptive Sparse PLS regression for continuous response
###    Short version for multiple call in cross-validation procedure
###
### Copyright 2014-10 Ghislain DURIF
###
### Adapted from R package "spls"
### Reference: Chun H and Keles S (2010)
### "Sparse partial least squares for simultaneous dimension reduction and variable selection",
### Journal of the Royal Statistical Society - Series B, Vol. 72, pp. 3--25.
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


spls.adapt.aux <- function(Xtrain, sXtrain, Ytrain, sYtrain, lambda.l1, ncomp, weight.mat, Xtest, sXtest, adapt=TRUE, meanXtrain, meanYtrain, sigmaXtrain, sigmaYtrain, center.X, center.Y, scale.X, scale.Y, weighted.center) {
	
	#####################################################################
	#### Initialisation
	#####################################################################
	
	Xtrain <- as.matrix(Xtrain)
	ntrain <- nrow(Xtrain) # nb observations
	p <- ncol(Xtrain) # nb covariates
	index.p <- c(1:p)
	Ytrain <- as.matrix(Ytrain)
	q <- ncol(Ytrain)
	
	sXtrain <- as.matrix(sXtrain)
	sYtrain <- as.matrix(sYtrain)
	
	Xtest <- as.matrix(Xtest)
	ntest <- nrow(Xtest)
	
	if(!is.null(weight.mat)) {
		V <- weight.mat
	} else {
		V <- diag(rep(1,ntrain), ncol=ntrain, nrow=ntrain)
	}
	
	
	
	#####################################################################
	#### Result objects
	#####################################################################
	betahat <- matrix(0, nrow=p, ncol=1)
	betamat <- list()
	X1 <- sXtrain
	Y1 <- sYtrain
	
	W <- matrix(data=NA, nrow=p, ncol=ncomp) # spls weight over each component
	T <- matrix(data=NA, nrow=ntrain, ncol=ncomp) # spls components
	P <- matrix(data=NA, nrow=ncomp, ncol=p) # regression of X over T
	Q <- matrix(data=NA, nrow=ncomp, ncol=q) # regression of Y over T
	
	
	#####################################################################
	#### Main iteration
	#####################################################################
	
	if ( is.null(colnames(sXtrain)) ) {
		Xnames <- index.p
	} else { 
		Xnames <- colnames(sXtrain)
	}
	
	new2As <- list()
	
	## SPLS
	for (k in 1:ncomp) {
		
		## define M
		M <- t(X1) %*% (V %*% Y1)
		
		#### soft threshold
		Mnorm1 <- median( abs(M) )
		
		M <- M / Mnorm1
		
		## adpative version
		if (adapt) {
			wi <- 1/abs(M)
			
			what <- ust.adapt(M, lambda.l1, wi)
			
		} else {
			## non adaptive version
			what <- ust(M, lambda.l1)
		}
		
		#### construct active set A
		A <- unique( index.p[ what!=0 | betahat[,1]!=0 ] )
		new2A <- index.p[ what!=0 & betahat[,1]==0 ]
		
		#### fit pls with selected predictors (meaning in A)
		X.A <- sXtrain[ , A, drop=FALSE ]           
		plsfit <- wpls( Xtrain=X.A, Ytrain=sYtrain, weight.mat=V, ncomp=min(k,length(A)), type="pls1", center.X=FALSE, scale.X=FALSE, center.Y=FALSE, scale.Y=FALSE, weighted.center=FALSE )
		
		
		#### output storage
		
		# weights
		w.k <- matrix(data=what, ncol=1)
		w.k <- w.k / sqrt(as.numeric(t(w.k) %*% w.k))
		W[,k] <- w.k
		
		# components on total observation space
		t.k <- (X1 %*% w.k) / as.numeric(t(w.k) %*% w.k)
		T[,k] <- t.k
		
		# regression of X over T
		p.k <- (t(X1) %*% (V %*% t.k)) / as.numeric(t(t.k) %*% (V %*% t.k))
		P[k,] <- t(p.k)
		
		# regression of Y over T
		q.k <- (t(Y1) %*% (V %*% t.k)) / as.numeric(t(t.k) %*% (V %*% t.k))
		Q[k,] <- t(q.k)
		
		
		## update
		
		Y1 <- sYtrain - plsfit$T %*% plsfit$Q
		X1 <- sXtrain
		X1[,A] <- sXtrain[,A] - plsfit$T %*% plsfit$P
		
		
		betahat <- matrix( 0, p, q )
		betahat[A,] <- matrix( plsfit$coeff, length(A), q )
		betamat[[k]] <- betahat # for cv.spls
		
		# variables that join the active set
		new2As[[k]] <- new2A
		
		
	}
	
	##### return objects
	
	hatYtest <- numeric(ntest)
	hatYtest.nc <- numeric(ntest)
	
	## predictions
	hatYtest <- sXtest %*% betahat
	
	#### betahat for non centered and non scaled data
	if((!scale.X) || (weighted.center)) { # if X non scaled, betahat don't have to be corrected regards sd.x
		sd.X <- rep(1, p)
	} else { # if X is scaled, it has to
		sd.X <- sigmaXtrain
	}
	if((!scale.Y) || (weighted.center)) {
		sd.Y <- 1
	} else {
		sd.Y <- sigmaYtrain
	}
	
	betahat.nc <- sd.Y * betahat / sd.X
	intercept <- meanYtrain - ( sd.Y * (drop( (meanXtrain / sd.X) %*% betahat)) )
	betahat.nc <- as.matrix(c(intercept, betahat.nc))
	
	#### non centered non scaled version of estimation and residuals
	hatYtest.nc <- cbind(rep(1,ntest),Xtest) %*% betahat.nc
	
	if ( !is.null(colnames(sXtrain)) ) {
		rownames(betahat) <- colnames(sXtrain)
	}
	
	#### return object
	result <- list( Xtrain=Xtrain, Ytrain=Ytrain, sXtrain=sXtrain, sYtrain=sYtrain, Xtest, sXtest,
				 betahat=betahat, betahat.nc=betahat.nc,
				 meanXtrain=meanXtrain, meanYtrain=meanYtrain, sigmaXtrain=sigmaXtrain, sigmaYtrain=sigmaYtrain,
				 hatYtest=hatYtest, hatYtest.nc=hatYtest.nc
				 )
	
	class(result) <- "spls.adapt.aux"
	return(result)
	
	
}
