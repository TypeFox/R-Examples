### wpls.R  (2014-10)
###
###    Weighted PLS regression
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



wpls <- function(Xtrain, Ytrain, ncomp, weight.mat=NULL, Xtest=NULL, type="pls1", center.X=TRUE, scale.X=FALSE, center.Y=TRUE, scale.Y=FALSE, weighted.center=FALSE) {
		
	#####################################################################
	#### Initialisation
	#####################################################################
	Xtrain <- as.matrix(Xtrain)
	ntrain <- nrow(Xtrain) # nb observations
	p <- ncol(Xtrain) # nb covariates
	index.p <- c(1:p)
	Ytrain <- as.matrix(Ytrain)
	q <- ncol(Ytrain)
	one <- matrix(1,nrow=1,ncol=ntrain)
	
	if(!is.null(Xtest)) {
		ntest <- nrow(Xtest)
	}
	
	# Tests on weighting matrix V
	if(!is.null(weight.mat)) { # weighting in scalar product (in observation space of dimension n)
		V <- as.matrix(weight.mat) 
		
		if ((!is.matrix(V)) || (!is.numeric(V))) {
			stop("Message from spls.adapt: V is not of valid type")}
		
		if ((ntrain != ncol(V)) || (ntrain != nrow(V))) {
			stop("Message from spls.adapt: wrong dimension for V, must be a square matrix of size the number of observations in Xtrain")
		}
	} else { # no weighting in scalar product
		V <- diag(rep(1, ntrain), nrow=ntrain, ncol=ntrain)
	}
		
		
	#####################################################################
	#### centering and scaling
	#####################################################################
	if (!weighted.center) {
		
		# Xtrain mean
		meanXtrain <- apply(Xtrain, 2, mean)
		
		# Xtrain sd
		sigmaXtrain <- apply(Xtrain, 2, sd)
		# test if predictors with null variance
		if ( any( sigmaXtrain < .Machine$double.eps )) {
			stop("Some of the columns of the predictor matrix have zero variance.")
		}
		
		# centering & eventually scaling X
		if(center.X && scale.X) {
			sXtrain <- scale( Xtrain, center=meanXtrain, scale=sigmaXtrain )
		} else if(center.X && !scale.X) {
			sXtrain <- scale( Xtrain, center=meanXtrain, scale=FALSE )
		} else {
			sXtrain <- Xtrain
		}
		
		# Y mean
		meanYtrain <- apply(Ytrain, 2, mean)
		
		# Y sd
		sigmaYtrain <- apply(Ytrain, 2, sd)
		# test if predictors with null variance
		if ( any( sigmaYtrain < .Machine$double.eps )) {
			stop("The response matrix has zero variance.")
		}
		# centering & eventually scaling Y
		if(center.Y && scale.Y) {
			sYtrain <- scale( Ytrain, center=meanYtrain, scale=sigmaYtrain )
		} else if(center.Y && !scale.Y) {
			sYtrain <- scale( Ytrain, center=meanYtrain, scale=FALSE )
		} else {
			sYtrain <- Ytrain
		}
		
		# Xtest
		if (!is.null(Xtest)) {
			
			## centering and scaling depend on Xtrain
			if(center.X && scale.X) {
				sXtest <- scale( Xtest, center=meanXtrain, scale=sigmaXtrain )
			} else if(center.X && !scale.X) {
				sXtest <- scale( Xtest, center=meanXtrain, scale=FALSE )
			} else {
				sXtest <- Xtest
			}
			
		}
		
	} else { # weighted scaling
		
		sumV <- sum(diag(V))
		
		# X mean
		meanXtrain <- matrix(diag(V), nrow=1) %*% Xtrain / sumV
		
		# X sd
		sigmaXtrain <- apply(Xtrain, 2, sd)
		# test if predictors with null variance
		if ( any( sigmaXtrain < .Machine$double.eps ) ) {
			stop("Some of the columns of the predictor matrix have zero variance.")
		}
		# centering & eventually scaling X
		sXtrain <- scale( Xtrain, center=meanXtrain, scale=FALSE )
		
		# Y mean
		meanYtrain <- matrix(diag(V), nrow=1) %*% Ytrain / sumV
		
		# Y sd
		sigmaYtrain <- apply(Ytrain, 2, sd)
		# test if predictors with null variance
		if ( any( sigmaYtrain < .Machine$double.eps ) ) {
			stop("The response matrix have zero variance.")
		}
		# centering & eventually scaling Y
		sYtrain <- scale( Ytrain, center=meanYtrain, scale=FALSE )
		
		# Xtest
		if (!is.null(Xtest)) {
			sXtest <- scale( Xtest, center=meanXtrain, scale=FALSE )
		}
		
	}
	
	
	#####################################################################
	#### Result objects
	#####################################################################
	
	X1 <- sXtrain
	Y1 <- sYtrain
	
	W <- matrix(data=NA, nrow=p, ncol=ncomp) # spls weight over each component
	T <- matrix(data=NA, nrow=ntrain, ncol=ncomp) # spls components
	P <- matrix(data=NA, nrow=ncomp, ncol=p) # regression of X over T
	Q <- matrix(data=NA, nrow=ncomp, ncol=q) # regression of Y over T
	
	
	#####################################################################
	#### Main iteration
	#####################################################################
	for ( k in 1:ncomp )
	{
		# direction vector
		
		w <- t(X1) %*% as.matrix( V %*% Y1 )
		w <- w / sqrt( as.numeric(t(w) %*% w) )
		W[,k] <- w
		
		# latent component
		t <- X1 %*% w / as.numeric(t(w) %*% w) 
		T[,k] <- t
		
		# standardize t?
		# t <- t - weighted.mean( t, diag(V) )
		# t <- t / sd(t)
		
		# coefficient 
		
		#coef.q <- t(t) %*% V %*% Y1 / drop( t(t) %*% V %*% t )
		coef.q <- (t(Y1) %*% (V %*% t)) / as.numeric( t(t) %*% (V %*% t) )
		Q[k,] <- t(coef.q)
		#coef.p <- t(t) %*% V %*% X1 / drop( t(t) %*% V %*% t )
		coef.p <- (t(X1) %*% (V %*% t)) / as.numeric( t(t) %*% (V %*% t) )
		P[k,] <- t(coef.p)
		
		# update
		if ( type=='pls1' )
		{
			Y1 <- Y1 - t %*% t(coef.q)
			X1 <- X1 - t %*% t(coef.p)
		}
		if ( type=='simpls' )
		{
			pj <- w
			pw <- pj %*% solve(t(pj) %*% pj) %*% t(pj)
			X1 <- X1 - X1 %*% pw
			#pw <- w %*% t(w) / drop( t(w) %*% w )
			#X1 <- X1 - X1 %*% pw
		}
	}
	
	
	## coeff B
	coeff <- W %*% Q
	betahat <- coeff
	
	#W.b <- W %*% solve( t(P) %*% W)
	#betahat <- W.b %*% t(Q)
	
	
	##### return objects
	
	## estim
	hatY <- sXtrain %*% betahat
	
	## residuals
	residuals <- sYtrain - hatY
	
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
	hatY.nc <- cbind(rep(1,ntrain),Xtrain) %*% betahat.nc
	residuals.nc <- Ytrain - hatY.nc
	
	## predictions
	if(!is.null(Xtest)) {
		hatYtest <- sXtest %*% betahat
		hatYtest.nc <- cbind(rep(1,ntest),Xtest) %*% betahat.nc
	} else {
		hatYtest <- NULL
		hatYtest.nc <- NULL
	}
	
	#### return object : list( W=W, T=T, Q=Q, P=P, coeff=coeff )
	result <- list( Xtrain=Xtrain, Ytrain=Ytrain, sXtrain=sXtrain, sYtrain=sYtrain,
				 betahat=betahat, betahat.nc=betahat.nc,
				 meanXtrain=meanXtrain, meanYtrain=meanYtrain, sigmaXtrain=sigmaXtrain, sigmaYtrain=sigmaYtrain,
				 X.score=T, X.loading=P, Y.loading=Q, X.weight=W, 
				 residuals=residuals, residuals.nc=residuals.nc,
				 hatY=hatY, hatY.nc=hatY.nc,
				 hatYtest=hatYtest, hatYtest.nc=hatYtest.nc,
				 ncomp=ncomp,
				 V=V, W=W, T=T, Q=Q, P=P, coeff=coeff)
	
	class(result) <- "wpls"
	return(result)
	
}