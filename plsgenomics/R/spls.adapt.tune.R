### spls.adapt.tune.R  (2014-10)
###
###    Tuning parameters (ncomp, lambda.l1) for adaptive sparse PLS regression for continuous response, by K-fold cross-validation
###
### Copyright 2014-10 Ghislain DURIF
###
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


spls.adapt.tune <- function(X, Y, lambda.l1.range, ncomp.range, weight.mat=NULL, adapt=TRUE, center.X=TRUE, center.Y=TRUE, scale.X=TRUE, scale.Y=TRUE, weighted.center=FALSE, return.grid=FALSE, ncores=1, nfolds=10) {
	
	#####################################################################
	#### Initialisation
	#####################################################################
	X <- as.matrix(X)
	n <- nrow(X) # nb observations
	p <- ncol(X) # nb covariates
	index.p <- c(1:p)
	Y <- as.matrix(Y)
	q = ncol(Y)
	one <- matrix(1,nrow=1,ncol=n)
	
	#####################################################################
	#### Tests on type input
	#####################################################################
	
	# On weighting matrix V
	if(!is.null(weight.mat)) { # weighting in scalar product (in observation space of dimension n)
		Vfull <- as.matrix(weight.mat) 
		
		if ((!is.matrix(Vfull)) || (!is.numeric(Vfull))) {
			stop("Message from spls.adapt: Vfull is not of valid type")}
		
		if ((n != ncol(Vfull)) || (n != nrow(Vfull))) {
			stop("Message from spls.adapt: wrong dimension for Vfull, must be a square matrix of size the number of observations in Xtrain")
		}
	} else { # no weighting in sclar product
		Vfull <- diag(rep(1,n), nrow=n, ncol=n)
	}
	
	# On weighted.center
	if ( (weighted.center) && (is.null(weight.mat))) {
		stop("Message from spls.adapt: if the centering is weighted, the weighting matrix V should be provided")
	}
	
	# ncores
	if ((!is.numeric(ncores)) || (round(ncores)-ncores!=0) || (ncores<1)) {
		stop("Message from rirls.spls.tune: ncores is not of valid type")
	}
	
	# nfolds
	if ((!is.numeric(nfolds)) || (round(nfolds)-nfolds!=0) || (nfolds<1)) {
		stop("Message from rirls.spls.tune: nfolds is not of valid type")
	}
	
	
	#####################################################################
	#### Cross-validation: computation on each fold over the entire grid
	#####################################################################
	
	## the train set is partitioned into nfolds part, each observation is assigned into a fold
	fold.obs <- sort(rep(1:nfolds, length.out = n))
	
	## hyper-parameter grid
	grid <- expand.grid(lambda.l1=lambda.l1.range, ncomp=ncomp.range, KEEP.OUT.ATTRS=FALSE)
	
	cv.grid.allfolds <- matrix( unlist( mclapply(1:nfolds, function(k) {
		
		
		#### train and test variable
		Xtrain <- subset(X, fold.obs != k)
		Ytrain <- subset(Y, fold.obs != k)
		
		ntrain <- nrow(Xtrain)
		
		Xtest <- subset(X, fold.obs == k)
		Ytest <- subset(Y, fold.obs == k)
		
		ntest <- nrow(Xtest)
		
		V <- Vfull[fold.obs == k, fold.obs == k]
		
		#####################################################################
		#### Tests on type input
		#####################################################################
		# On Xtrain
		if ((!is.matrix(Xtrain)) || (!is.numeric(Xtrain))) {
			stop("Message from spls.adapt: Xtrain is not of valid type")
		}
		
		if (p==1) {
			stop("Message from spls.adapt: p=1 is not valid")}
		
		# On Xtest if necessary	
		if (is.vector(Xtest)==TRUE) {
			Xtest <- matrix(Xtest,nrow=1)
		}
		
		Xtest <- as.matrix(Xtest)
		ntest <- nrow(Xtest) 
		
		if ((!is.matrix(Xtest)) || (!is.numeric(Xtest))) {
			stop("Message from spls.adapt: Xtest is not of valid type")}
		
		if (p != ncol(Xtest)) {
			stop("Message from spls.adapt: columns of Xtest and columns of Xtrain must be equal")
		}
		
		# On Ytrain
		if ((!is.matrix(Ytrain)) || (!is.numeric(Ytrain))) {
			stop("Message from spls.adapt: Ytrain is not of valid type")
		}
		
		if (q != 1) {
			stop("Message from spls.adapt: Ytrain must be univariate")
		}
		
		if (nrow(Ytrain)!=ntrain) {
			stop("Message from spls.adapt: the number of observations in Ytrain is not equal to the Xtrain row number")
		}
		
		# On hyper parameter: lambda.ridge, lambda.l1
		if ( (sum(!is.numeric(lambda.l1.range))) || (sum(lambda.l1.range<0)) || (sum(lambda.l1.range>1)) ) {
			stop("Message from spls.adapt: lambda is not of valid type")
		}
		
		# ncomp type
		if ( (sum(!is.numeric(ncomp.range))) || (sum(round(ncomp.range)-ncomp.range!=0)) || (sum(ncomp.range<1)) || (sum(ncomp.range>p)) ) {
			stop("Message from spls.adapt: ncomp is not of valid type")
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
				sXtrain <- scale( Xtrain, center=meanXtrain, scale=sigmaXtrain)
			} else if(center.X && !scale.X) {
				sXtrain <- scale( Xtrain, center=meanXtrain, scale=FALSE)
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
				sYtrain <- scale( Ytrain, center=meanYtrain, scale=sigmaYtrain)
			} else if(center.Y && !scale.Y) {
				sYtrain <- scale( Ytrain, center=meanYtrain, scale=FALSE)
			} else {
				sYtrain <- Ytrain
			}
			
			if(center.Y && scale.Y) {
				sYtest <- scale( Ytest, center=meanYtrain, scale=sigmaYtrain)
			} else if(center.Y && !scale.Y) {
				sYtest <- scale( Ytest, center=meanYtrain, scale=FALSE)
			} else {
				sYtest <- Ytest
			}
			
			# Xtest	
			## centering and scaling depend on Xtest
			if(center.X && scale.X) {
				sXtest <- scale( Xtest, center=meanXtrain, scale=sigmaXtrain )
			} else if(center.X && !scale.X) {
				sXtest <- scale( Xtest, center=meanXtrain, scale=FALSE )
			} else {
				sXtest <- Xtest
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
			sYtest <- scale( Ytest, center=meanYtrain, scale=FALSE )
			
			# Xtest
			sXtest <- scale( Xtest, center=meanXtrain, scale=FALSE )
			
		}
		
		
		#####################################################################
		#### Computation over the entire grid
		#####################################################################
		
		### sur chaque fold, on calcule pour tout lambda.ridge.range
		cv.grid.byfold <- matrix( sapply( split(grid, f=row.names(grid)), function(grid.line) {
			
			model <- tryCatch( spls.adapt.aux(Xtrain=Xtrain, sXtrain=sXtrain, Ytrain=Ytrain, sYtrain=sYtrain, lambda.l1=grid.line$lambda.l1, ncomp=grid.line$ncomp, weight.mat=weight.mat, Xtest=Xtest, sXtest=sXtest, adapt=adapt, meanXtrain=meanXtrain, meanYtrain=meanYtrain, sigmaXtrain=sigmaXtrain, sigmaYtrain=sigmaYtrain, center.X=center.X, center.Y=center.Y, scale.X=scale.X, scale.Y=scale.Y, weighted.center=weighted.center), error = function(e) { warnings("Message from spls.adapt.tune: error when fitting a model in crossvalidation"); return(NULL);} )
			
			## resutls
			res <- numeric(4)
			
			if(!is.null(model)) {
				res <- c(grid.line$lambda.l1, grid.line$ncomp, k, sum((model$hatYtest - sYtest)^2) / ntest)
			} else {
				res <- c(grid.line$lambda.l1, grid.line$ncomp, k, NA)
			}
			
			return(res)
			
		}), ncol=4, byrow=TRUE)
		
		return( t(cv.grid.byfold) )
		
		
	}, mc.cores = ncores, mc.silent=TRUE)), ncol=4, byrow=TRUE)
	
	cv.grid.allfolds <- data.frame(cv.grid.allfolds)
	colnames(cv.grid.allfolds) <- c("lambda.l1", "ncomp", "nfold", "error")
	
	#####################################################################
	#### Find the optimal point in the grid
	#####################################################################
	
	
	## compute the mean error over the folds for each point of the grid
	cv.grid <- data.frame( as.matrix( with( cv.grid.allfolds, aggregate(error, list(lambda.l1, ncomp), function(x) c(mean(x, na.rm=TRUE), sd(x, na.rm=TRUE)) ))))
	colnames(cv.grid) <- c("lambda.l1", "ncomp", "error", "error.sd")
	
	
	##### return
	if(return.grid) {
		
		return( list(lambda.l1.opt = cv.grid$lambda.l1[which.min(cv.grid$error)], ncomp.opt = cv.grid$ncomp[which.min(cv.grid$error)], cv.grid=cv.grid) )
		
	} else {
		
		return( list(lambda.l1.opt = cv.grid$lambda.l1[which.min(cv.grid$error)], ncomp.opt = cv.grid$ncomp[which.min(cv.grid$error)], cv.grid=NULL) )
		
	}
	
	
}
	
