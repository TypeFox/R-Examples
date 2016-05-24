### rirls.spls.tune.R  (2014-10)
###
###    Tuning parameters (ncomp, lambda.l1, lambda.ridge) for Ridge Iteratively Reweighted Least Squares followed by Adaptive Sparse PLS regression for binary response, by K-fold cross-validation
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


rirls.spls.tune <- function(X, Y, lambda.ridge.range, lambda.l1.range, ncomp.range, adapt=TRUE, maxIter=100, svd.decompose=TRUE, return.grid=FALSE, ncores=1, nfolds=10) {
	
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
	
	# maxIter
	if ((!is.numeric(maxIter)) || (round(maxIter)-maxIter!=0) || (maxIter<1)) {
		stop("Message from rirls.spls.tune: maxIter is not of valid type")
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
	fold.obs = sort(rep(1:nfolds, length.out = n))
	
	## hyper-parameter grid
	grid = expand.grid(lambda.ridge=lambda.ridge.range, lambda.l1=lambda.l1.range, ncomp=ncomp.range, KEEP.OUT.ATTRS=FALSE)
	
	cv.grid.allfolds = matrix( unlist( mclapply(1:nfolds, function(k) {
		
		
		#### train and test variable
		Xtrain <- subset(X, fold.obs != k)
		Ytrain <- subset(Y, fold.obs != k)
		
		ntrain <- nrow(Xtrain)
		
		Xtest <- subset(X, fold.obs == k)
		Ytest <- subset(Y, fold.obs == k)
		
		ntest <- nrow(Xtest)
		
		#####################################################################
		#### Tests on type input
		#####################################################################
		# On Xtrain
		if ((!is.matrix(Xtrain)) || (!is.numeric(Xtrain))) {
			stop("Message from rirls.spls.tune: Xtrain is not of valid type")
		}
		
		if (p==1) {
			stop("Message from rirls.spls.tune: p=1 is not valid")}
		
		# On Xtest
		if (is.vector(Xtest)==TRUE) {
			Xtest <- matrix(Xtest,nrow=1)
		}
		
		Xtest <- as.matrix(Xtest)
		ntest <- nrow(Xtest) 
		
		if ((!is.matrix(Xtest)) || (!is.numeric(Xtest))) {
			stop("Message from rirls.spls.tune: Xtest is not of valid type")}
		
		if (p != ncol(Xtest)) {
			stop("Message from rirls.spls.tune: columns of Xtest and columns of Xtrain must be equal")
		}	
		
		# On Ytrain
		if ((!is.matrix(Ytrain)) || (!is.numeric(Ytrain))) {
			stop("Message from rirls.spls.tune: Ytrain is not of valid type")
		}
		
		if (q != 1) {
			stop("Message from rirls.spls.tune: Ytrain must be univariate")
		}
		
		if (nrow(Ytrain)!=ntrain) {
			stop("Message from rirls.spls.tune: the number of observations in Ytrain is not equal to the Xtrain row number")
		}
		
		# On Ytrain value
		if (sum(is.na(Ytrain))!=0) {
			stop("Message from rirls.spls.tune: NA values in Ytrain")
		}
		
		if (sum(!(Ytrain %in% c(0,1)))!=0) {
			stop("Message from rirls.spls.tune: Ytrain is not of valid type")
		}
		
		if (sum(as.numeric(table(Ytrain))==0)!=0) {
			stop("Message from rirls.spls.tune: there are empty classes")
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
				stop("Message from rirls.spls.tune: the procedure stops because number of predictor variables with no null variance is less than 1.")
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
		
		# center and scale Xtest	
		meanXtest <- apply(Xtest,2,mean)
		sigma2test <- apply(Xtest,2,var)
		
		sXtest <- scale(Xtest, center=meanXtrain, scale=sqrt(sigma2train))
		
		sXtest.nosvd <- sXtest # keep in memory if svd decomposition
		
		# if svd decomposition
		if ((p > ntrain) && (svd.decompose)) {
			sXtest <- sXtest%*%V
		}
		
		
		
		#####################################################################
		#### Computation over the entire grid
		#####################################################################
		
		### sur chaque fold, on calcule pour tout lambda.ridge.range
		cv.grid.byfold <- matrix( sapply( split(grid, f=row.names(grid)), function(grid.line) {
			
			# On hyper parameter: lambda.ridge, lambda.l1
			if ((!is.numeric(grid.line$lambda.ridge)) || (grid.line$lambda.ridge<0) || (!is.numeric(grid.line$lambda.l1)) || (grid.line$lambda.l1<0)) {
				stop("Message from rirls.spls.tune: lambda is not of valid type")
			}
			
			# ncomp type
			if ((!is.numeric(grid.line$ncomp)) || (round(grid.line$ncomp)-grid.line$ncomp!=0) || (grid.line$ncomp<0) || (grid.line$ncomp>p)) {
				stop("Message from rirls.spls.tune: ncomp is not of valid type")
			}
			
			model <- tryCatch( rirls.spls.aux(sXtrain=sXtrain, sXtrain.nosvd=sXtrain.nosvd, Ytrain=Ytrain, lambda.ridge=grid.line$lambda.ridge, lambda.l1=grid.line$lambda.l1, ncomp=grid.line$ncomp, sXtest=sXtest, sXtest.nosvd=sXtest.nosvd, adapt=adapt, maxIter=maxIter, svd.decompose=svd.decompose, meanXtrain=meanXtrain, sigma2train=sigma2train), error = function(e) { warnings("Message from rirls.spls.tune: error when fitting a model in crossvalidation"); return(NULL);} )
			
			## resutls
			res = numeric(6)
			
			if(!is.null(model)) {
				res = c(grid.line$lambda.ridge, grid.line$lambda.l1, grid.line$ncomp, k, model$converged, sum(model$hatYtest != Ytest) / ntest)
			} else {
				res = c(grid.line$lambda.ridge, grid.line$lambda.l1, grid.line$ncomp, k, NA, NA)
			}
			
			return(res)
			
		}), ncol=6, byrow=TRUE)
		
		return( as.vector(t(cv.grid.byfold)) )
		
		
	}, mc.cores=ncores, mc.silent=TRUE)), ncol=6, byrow=TRUE)
	
	cv.grid.allfolds = data.frame(cv.grid.allfolds)
	colnames(cv.grid.allfolds) = c("lambda.ridge", "lambda.l1", "ncomp", "nfold", "converged", "error")
	
	#####################################################################
	#### Find the optimal point in the grid
	#####################################################################
	
	
	## compute the mean error over the folds for each point of the grid
	cv.grid.error <- data.frame( as.matrix( with( cv.grid.allfolds, aggregate(error, list(lambda.ridge, lambda.l1, ncomp), function(x) c(mean(x, na.rm=TRUE), sd(x, na.rm=TRUE)) ))))
	colnames(cv.grid.error) <- c("lambda.ridge", "lambda.l1", "ncomp", "error", "error.sd")
	
	## compute the percentage of convergence over the folds for each point of the grid
	cv.grid.conv <- data.frame( as.matrix( with( cv.grid.allfolds, aggregate(converged, list(lambda.ridge, lambda.l1, ncomp), mean, na.rm=TRUE))))
	colnames(cv.grid.conv) <- c("lambda.ridge", "lambda.l1", "ncomp", "converged")
	
	## % of convergence over all cross-validation process
	conv.per <- mean(cv.grid.conv$converged)
	
	## merge the two tables
	cv.grid <- merge(cv.grid.error, cv.grid.conv, by = c("lambda.ridge", "lambda.l1", "ncomp"))
	
	##### return
	if(return.grid) {
		
		return( list(lambda.ridge.opt = cv.grid$lambda.ridge[which.min(cv.grid$error)], lambda.l1.opt = cv.grid$lambda.l1[which.min(cv.grid$error)], ncomp.opt = cv.grid$ncomp[which.min(cv.grid$error)], conv.per=conv.per, cv.grid=cv.grid) )
		
	} else {
		
		return( list(lambda.ridge.opt = cv.grid$lambda.ridge[which.min(cv.grid$error)], lambda.l1.opt = cv.grid$lambda.l1[which.min(cv.grid$error)], ncomp.opt = cv.grid$ncomp[which.min(cv.grid$error)], conv.per=conv.per, cv.grid=NULL) )
		
	}
	
	
	
	
	
	
	
}
	
