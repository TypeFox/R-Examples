###
# gensemble
###

setClass("gensemble",
	representation(
		abm="AbstractModel",
		dclass="logical",
		nlev="numeric",
		ylevels="vector",
		mods="list",
		nmods="numeric",
		bagmat="matrix",
		oobpred="ANY", #a vector or a matrix
		oobpredmat="matrix",
		accmat="matrix",
		test_oobpred="ANY",
		test_oobpredmat="matrix",
		test_accmat="matrix"
	),
)

###
# print()
# pretty minimalist, needs more detail
###
print.gensemble <- function(x, ...)
{ 
	lb <- x
	if (!is.null(lb@abm)) {
		cat(paste('Model Type:', lb@abm@model), '\n')
	} else {
		cat('gensemble: no model defined\n')
		return
	}
	if (lb@nmods > 0) {
		cat(sprintf('%d models built', length(lb@mods)), '\n')
	}
}

###
# returns the proportions of each class in a given vector
###
lb.class.prop <- function(rowX, labs)
{
	p <- c()
	for (lev in labs) {
		p <- c(p, length(rowX[rowX==lev])/length(rowX))
	}
	return(t(p))	
}

###
# the predict() function
###
predict.gensemble <- function(object, X, type=c("prob", "class"), method=c("prob", "vote"), return.all=F, ...)
{
	LB <- object
	abm <- LB@abm
	mods <- LB@mods
	dclass <- LB@dclass
	
	if (dclass) {
		#classification
		nlev <- LB@nlev
		labs <- LB@ylevels
		method <- match.arg(method, c("prob", "vote"))
		type <- match.arg(type, c("prob", "class"))

		pred <- matrix(0, nrow(X), nlev)
		if (method == 'vote')
			predmat <- matrix(0, nrow(X), length(mods))	
	} else {
		#regression
		pred <- rep(0, nrow(X))
		if (return.all)
			predmat <- matrix(0, nrow(X), length(mods))
	}

	for (i in 1:length(mods)) {
		newpred <- ab.predict(abm, mods[[i]], X)
		if (dclass && method == 'vote') {
			predmat[,i] <- apply(newpred, 1, class.or.na, labels=labs)
		} else {
			pred <- pred + newpred
			if (return.all)
				predmat[,i] <- newpred
		}
	}
	
	if (dclass && method == 'vote') {
		#get the prop of votes for that level
		pred <- apply(predmat, 1, lb.class.prop, labs=labs)
		pred <- t(pred)
	} else {
		pred <- pred / length(mods)
	}

	if (dclass && type == 'class') {
		return(as.factor(apply(pred, 1, class.or.na, labs)))
	} else {
		if (!dclass && return.all)
			return(list(pred, predmat))
		return(pred)
	}	
}


###
# linbag(model, X, Y, ...)
# todo: check weights
###
gensemble <- function(abm, X, Y, sampsize=NULL, sampsize_prop=FALSE, nmods=100, perturb_val=0.1, Xtest=NULL, Ytest=NULL, do.trace=TRUE, stepsize=10)
{
	
	#if we were passed a test set, do some sanity checking
	if (!is.null(Xtest)) {
		if (nrow(Xtest) != length(Ytest)) 
			stop("test data set size mismatch")
		testset <- TRUE
	} else {
		testset <- FALSE
	}
	
	#if we have a factor in Y, we are doing classification and keep a different set of tracking metrics
	if (is.factor(Y)) {
		dclass <- TRUE
		nlev <- nlevels(Y)
		#per iteration summary function
		#this is tied to what mk_oobpred_storage returns as well
		sfunc <- cls_iter_summary
	
	} else {
		#doing regression
		dclass <- FALSE
		nlev <- 1
		sfunc <- reg_iter_summary	
	}

	#bagging matrix, keeps track of which obs used in which iteration
	bagmat <- matrix(FALSE, nrow(X), nmods)
	#list of models generated
	modlist <- list(NULL)
	
	#Out of bag prediction matrix. Used to store the predictions of the model from each iteration
	#used to build the out of bag predictions
	tmp <- mk_oobpred_storage(Y, nmods)
	oobpredmat <- tmp[[1]]
	#our current best prediction set
	oobpred <- tmp[[2]]
	#tracking the accuracy per iteration
	accmat <- tmp[[3]]
	#if we were passed a test set we will keep track of how it goes as well
	if (testset) {
		tmp <- mk_oobpred_storage(Ytest, nmods)
		test_oobpredmat <- tmp[[1]]
		test_oobpred <- tmp[[2]]
		test_accmat <- tmp[[3]]
	}
	
	#generate sample sizes for each iteration, this will handle null sampsize
	#and try to come up with something sensible for classification and regression
	sizes <- mksampsize(Y, sampsize, sampsize_prop)
	
	for (i in 1:nmods) {

		#generate the indexes to use for this iteration based on sample size(s)
		bagsamp <- bagsamp_idx(Y, sizes, rep=T)
		bagmat[bagsamp, i] <- TRUE
		Yout <- Y[bagsamp]
		Yout <- perturb(Yout, perturb_val)
		
		sink("modelouts.txt")
		modlist[[i]] <- ab.model(abm, X[bagsamp,], Yout)
		sink()
		
		#which cols in the oob prediction matrix to update
		starti <- (i - 1) * nlev + 1
		endi <- starti + nlev - 1
		
		#set the col(s) to 0
		oobpredmat[,starti:endi] <- rep(NA, length(Y))
		#run the new model on the out of bag data
		oobpredmat[-bagsamp, starti:endi] <- ab.predict(abm, modlist[[i]], X[-bagsamp,])
		#see how the model does on the test set
		if (testset)
			test_oobpredmat[,starti:endi] <- ab.predict(abm, modlist[[i]], Xtest)
		
		#update oobpred with the mean values for all model results stored in oobpredmat
		if (dclass) {
			#get the columns of the prediction matrix for a given level
			#aggregate the predictions in oobpredmat using the agg_oobpred function 
			#store this in oobpred
			for (lev in levels(Y)) {
				oobpred[,levcols(oobpred,lev)] <- apply(oobpredmat[,levcols(oobpredmat, lev)], 1, agg_oobpred)
				if (testset)
					test_oobpred[,levcols(test_oobpred, lev)] <- apply(test_oobpredmat[,levcols(test_oobpredmat, lev)], 1, agg_oobpred)
			}
		} else {
			oobpred <- apply(oobpredmat, 1, agg_oobpred)
			if (testset)
				test_oobpred <- apply(test_oobpredmat, 1, agg_oobpred)
		}
#		cat(oobpredmat[,starti:endi], '\n')
#		cat(oobpredmat, '\n')
#		cat(oobpred, '\n')

		#book keeping time
		accmat[i, ] <- sfunc(oobpred, Y)
		if (testset) 
			test_accmat[i, ] <- sfunc(test_oobpred, Ytest)	
			
		#print status
		if(do.trace && i %% stepsize == 0) {
			dimnm <- dimnames(accmat)
			if (i == stepsize) {
				cat(formatC(dimnm[[2]], format="s", width=10))
				if (testset) {
					testnames <- paste("test", dimnm[[2]], sep=".")	
					testwidth <- max(nchar(testnames)) 				
					cat(' ', formatC(testnames, format="s", width=testwidth))
				}
				cat('\n')
			}
				
			pvec <- formatC(as.vector(accmat[i,]), width=10, digits=4, format="f")
			if (testset) 
				pvec <- c(pvec, formatC(as.vector(test_accmat[i,]), width=testwidth, digits=4, format="f"))
			cat(pvec, "\n")
		}
	}

	#create a new linboost object	
	lb <- new("gensemble", abm=abm)
	if (dclass) {
		lb@dclass <- TRUE
		lb@nlev <- nlev
		lb@ylevels <- levels(Y)
	} else {
		lb@dclass <- FALSE
	}
	
	lb@mods <- modlist
	lb@nmods <- length(modlist)
	lb@bagmat <- bagmat
	lb@oobpred <- oobpred
	lb@oobpredmat <- oobpredmat
	lb@accmat <- accmat
	if (testset) {
		lb@test_oobpred <- test_oobpred
		lb@test_oobpredmat <- test_oobpredmat
		lb@test_accmat <- test_accmat
	}

	return(lb)	

###
# old return, delete this before release
###
#	if (testset) {
#		LB <- list(modlist, bagmat, oobpred, oobpredmat, accmat, test_oobpred, test_oobpredmat, test_accmat)
#	} else {
#		LB <- list(modlist, bagmat, oobpred, oobpredmat, accmat)
#	}
#	return(LB)

}

###
# the summary information collected for each iteration of regression
# changes to this should also change mk_oobpred_storage
###
reg_iter_summary <- function(yhat, y)
{
	#for regression, we want the mse, variance, mse scaled by variance
	vr <- var(yhat, na.rm=T)
	mse <- calc_mse(yhat, y)
	scaled_mse <- mse / vr
	psuedo_r2 <- cor(y, yhat, use="complete.obs")
	na <- sum(is.na(yhat))/length(y)
	c(vr, mse, scaled_mse, psuedo_r2, na)
}

###
# the summary information collected for each iteration of classification
# changes to this should also change mk_oobpred_storage
###
cls_iter_summary <- function(yhat, y)
{
	#XXX TODO: fix these names, should be class.hit etc
	s <- vector("numeric", length=0)
	labs <- apply(yhat, 1, class.or.na, labels=levels(y))
	#per class accuracy
	for (lev in levels(y)) {
		s <- c(s,pred_error(y, labs, lev)["hit"])
	}
	#proportion that have not yet had classification attempted
	s <- c(s, sum(is.na(labs)) / length(labs))
	#total accuracy
	s <- c(s, pred_error(y, labs)["hit"])
	return(s)
}

###
# creates storage for the internal operations of linboost
# oobpredmat keeps track of what was predicted for each iteration
# oobpred is the running "best" prediction updated every iteration
# accmat stores information about accuracy of each iteration
#
# the specific columns differ for regression and classfication, so
# are generated by reg_iter_summary and cls_iter_summary respectively
###
mk_oobpred_storage <- function(Y, nmods) 
{
	if (!is.factor(Y)) {
		oobpredmat <- matrix(NA, length(Y), nmods)
		oobpred <- rep(0, length(Y))
		accmat <- matrix(0, nmods, 5)
		colnames(accmat) <- c("var", "mse", "scaled.mse", "pseudo R^2", "obs.na")		
	} else {
		#track accuracy over iterations, +2 to track unlabelled and overall accuracy		
		nlev <- length(levels(Y))
		#for classification, we store a column for each level defined
		oobpredmat <- matrix(NA, length(Y), nlev * nmods)
		colnames(oobpredmat) <- rep(levels(Y), nmods)
		#the out of bag prediction
		oobpred <- matrix(0, length(Y), nlev)
		colnames(oobpred) <- levels(Y)
		accmat <- matrix(0, nmods, nlev+2)
		colnames(accmat) <- c(levels(Y), "na", "total")	
	}
	return(list(oobpredmat, oobpred, accmat))
}

###
# come up with some meaningful intepretation of the sample size 
# provided to linboost. this is used when creating the bagged sample.
# if proporption is TRUE, all sample sizes provided are assumed to be
# proportions and should range between 0 and 1
#
# for regression if the sample size is null, it is assumed to be the length
# of the input data.
#
# if the given Y value is not a factor, regression is assumed, both here and 
# in linboost
#
# for classification, we want a list keyed by each level with the value being
# the sample size for that class.
#
# if sampsize is null, it will build a list of the levels with the count for
# each level
#
# if it is a list, and proportion is false, it will check the provided levels
# match what is available in Y, and return the sample values passed in with no
# modification
#
# if it is a list, and proportion is true, it will set the value for each level 
# to be the proportion specifed. So if Level A has n 100 and sampsize[["A"]] is 0.7
# the returned list will be sampesize[["A"]] = 70
#
# if it is a vector, it must either have 1 element or an element for each level. 
# 
# if there is 1 element, it is assumed to be the sample size for each level. 
#
# If proportion is TRUE, the sample size will be generated based on the count 
# of each level in Y. E.g. if Y has 100 of level A and 50 of level B, sampsize
# is 0.7 and proportion is TRUE, the output would be list(A=70, B=35)
#
# If there is 1 element and proportion is false, it will be the sample size for
# each level.
#
# If it is a vector of length greater than 1, it is assumed each item corresponds
# the the order of levels in Y as returned by levels(). If you want to guarantee
# the level names match the right size, pass a list of the level names and the 
# sizes you want. (see above for details on this)
#
# note the actual sampling is done by bagsamp_idx, which requires a valid
# sampsize argument. use this function to generate the valid sampsize argument
mksampsize <- function(Y, sampsize=NULL, proportion=FALSE) 
{
	#for regression it is either specified or we use the length of the input
	if (!is.factor(Y)) {
		if (is.null(sampsize)) {
			sz <-round(length(Y) * 0.8)
		} else {
			if (proportion) {
				sz <- round(length(Y) * sampsize)
			} else {
				sz <- sampsize
			}
		}
		return(sz)
	}

	#for classification, we use a list with a sample size for each class
	
	sz <- list()
	
	if (is.null(sampsize)) {
		#no sample size given so use the counts per level
		for (lev in levels(Y)) {
			sz[[lev]] <- round(length(which(Y == lev) * 0.8))
		}
	} else if (is.list(sampsize)) {
		#a list of sample sizes per class
		#check the names match
		for (nm in names(sampsize)) {
			if (! nm %in% levels(Y)) {
				stop("unknown level name in sampsize")
			}
		}
		
		if (proportion) {
			for (lev in levels(Y)) {
				sz[[lev]] <- round(length(which(Y == lev)) * sampsize[lev])
			}
		} else {
			sz <- sampsize
		}
	} else if (is.vector(sampsize)) {
		#a vector of sample sizes given, assume it matches the order of levels(Y)
		#if we were passed a single value, assume it is the same value for each level
		if (length(sampsize) == 1) {
			sampsize <- rep(sampsize, length(levels(Y)))
		}
		
		if (length(sampsize) != length(levels(Y))) {
			stop("sample size vector does not match number of levels")
		}
		for (i in 1:length(sampsize)) {
			lev <- levels(Y)[i]
			if (proportion) {
				sz[[lev]] <- round(length(which(Y == lev)) * sampsize[i])
			} else { 
				sz[[lev]] <- sampsize[i]
			}
		}
	} else {
		stop("dont know what to do with sampsize")
	}
	return(sz)
}

###
# sample the given Y values based on the sampsize argument
# 
# for regression, sampsize should be the length of the sample required.
#
# for classification, it should be a list keyed by each class in Y with 
# a corresponding value that is the desired sample size for that class.
#
# see mksampsize for details the sampsize argument
#
# XXX todo weights.
###
bagsamp_idx <- function(Y, sampsize, weights=NULL, rep=T)
{
	#assumes sampsize has been set up correctly via mksampsize
	if (!is.factor(Y)) {
		return(sample(1:length(Y), sampsize, replace=rep, prob=weights))
	}
	#need to check how weights work, probably a weight per class
	bagsamp <- c()
	for (lev in levels(Y)) {
		idx <- which(Y == lev)
		sz <- sampsize[[lev]]
		bagsamp <- c(bagsamp, sample(idx, sz, replace=rep, prob=weights))
	}	
	return(bagsamp)
}

###
# perturb the inputs
###
perturb <- function(y, prop=0) {
	if (prop == 0L){
		#noop
		return(y)
	} else if (prop > 1) {
		stop("perturb proportion should be between 0 and 1")
	}
	
	#how many items in y to perturb
	plen <- round(length(y) * prop)
	#inidices to perturb
	pvec <- sample(1:length(y), plen)
	#do the peturb
	y[pvec] <- y[pvec][order(rnorm(plen))]
	return(y)
}

###
# calculate the prediction error for classification
# yt is the true y values, yhat are the predicted values
# pass level to check for a given level accuracy only.
##
# xxx swap yt and yhat args to be consistent w/other funcs
###
pred_error <- function(yt, yhat, level=NULL) {
	
	if (is.null(level)) {
		#look at the whole set
		y.idx <- 1:length(yt)
	} else {
		#look at a specific level
		y.idx <- which(yt == level)
	}
	
	if (length(y.idx) == 0) {
		return(NA)
	}
	
	# hit = proportion of correct classification
	hit <- sum(na.omit(yhat[y.idx] == yt[y.idx])) / length(na.omit(yhat[y.idx]))
	# miss = proportion of incorrect classification
	miss <- sum(na.omit(yhat[y.idx] != yt[y.idx])) / length(na.omit(yhat[y.idx]))
	#stopifnot(miss == 1 - hit)
	c(hit=hit, miss=miss)
}

###
# calculate the mean squared error
###
calc_mse <- function(yhat, y) {
	return(sum((yhat - y)^2, na.rm=T)/ length(y))
}

###
# return the most likely class based on the x values
# returns NA if NaN's exist in X
###
class.or.na <- function(X, labels) {
	if (is.nan(sum(X))) {
		return(NA)
	}
	return(labels[which.max(X)])
}

###
#take the mean of all values, skipping NA's
###
agg_oobpred <- function(X) {
	mean(X, na.rm=T)
}

###
# return the column indexes that match the given level. used when
# aggregating prediction summaries.
###
levcols <- function(mat, level) {
	which(colnames(mat) == level)
}
