### EMtree-functions.r
# This file contains helper functions for the 
# EMtree.r code.


### Function: tree
# This is the generic (S3-type) of function that returns the tree associated
# with an object (such as a RE-EM tree or a FE-EM tree).
tree <- function(object,...){
    UseMethod("tree")
}


	
### Function: REEMtree
# This program runs the main EM/tree algorithm of
# Simonoff & Sela.
# 6/5/08 (RS): To deal with subsets, we now subset the entire dataset,
# use that for fitting, and then deal with any random effects that
# were not estimated by setting them to 0 and printing a warning.
# 6/14/08 (RS): Default method is now REML, because ML doesn't always converge.
# 7/31/08 (RS): Added in a few lines to deal with the case when the tree doesn't split,
# since LME can't run with a single level of contrasts.
# 6/24/09 (RS): No longer need to identify the grouping variable separately (but there
# can still be only one grouping variable.
# 6/25/09 (RS): Replaced SubtractEffects and just use the difference between the residuals
# at the different levels to compute the random effects.
# Inputs:
#	formula - a formula (as in the usual LM/rpart calls)
#	data - data (as in the usual LM/rpart calls)
#	random - the formula describing the random effects (as in LME)
#	groups - a vector containing the group identifier for each observation
#		(in the same order as the observations) - this is identical to
#		the variable after | in random, but we can't extract that directly
#		for some reason)
#	subset - subset of the data to use for fitting
#	initialRandomEffects [0] - a vector of initial values for random effects
#	ErrorTolerance [0.001]
#	MaxIterations [1000]
# 	likelihoodCheck [TRUE] - should the likelihood be used to check for
#		convergence?  (if not, the random effects are checked instead)
#	verbose [FALSE] - print intermediate trees
#	### These options pertain to the RPART part of estimation
#	tree.control [rpart.control] - controls to be passed through to rpart
#	cv [TRUE] - Should cross-validation be used?
#	cpmin [0.0001] - complexity parameter used in building a tree before cross-validation
#	cpcv [0.01] - complexity used for pruning in a cross-validated tree
#	no.SE [0] - number of standard errors used in pruning (0 if unused)
#	### These options pertain to the LME part of estimation
#	lme.control [lmeControl(returnObject=TRUE)] - controls to be passed through to LME
#	method ["REML"] - "ML" or "REML", depending on whether the effects should be estimated
#		with maximum likelihood or restricted maximum likelihood
#	correlation [NULL] - an option CorStruct object describing the within-group correlation structure
REEMtree <- function(formula, data, random, subset=NULL, initialRandomEffects=rep(0,TotalObs), 
		ErrorTolerance=0.001, MaxIterations=1000, verbose=FALSE, tree.control=rpart.control(), 
		cv=TRUE, cpmin = 0.001, no.SE =1,
		lme.control=lmeControl(returnObject=TRUE), method="REML", correlation=NULL){
    TotalObs <- dim(data)[1]

    originaldata <- data

    # Subset the data if necessary
    if(identical(subset, NULL)){ 
	subs <- rep(TRUE, dim(data)[1])
    } else {
	subs <- subset
    }

    # Parse formula
    Predictors <- paste(attr(terms(formula),"term.labels"),collapse="+")
    TargetName <- formula[[2]]
    # Remove the name of the data frame if necessary
    if(length(TargetName)>1) 
	TargetName <-TargetName[3]
    if(verbose) print(paste("Target variable: ", TargetName))

    Target <- data[,toString(TargetName)]


    # Condition that indicates the loop has not converged or 
    # run out of iterations
    ContinueCondition <- TRUE

    iterations <- 0

    # Get initial values
    AdjustedTarget <- Target - initialRandomEffects
    oldlik <- -Inf

    # Make a new data frame to include all the new variables
    newdata <- data
    newdata[, "SubsetVector"] <- subs

    while(ContinueCondition){

	# Current values of variables
  	newdata[,"AdjustedTarget"] <- AdjustedTarget
	iterations <- iterations+1

	# Compute current tree
         if (cv) {
             tree1 <- rpart(formula(paste(c("AdjustedTarget", Predictors),
                 collapse = "~")), data = newdata, subset = subs,
                 method = "anova", control = rpart.control(cp=cpmin))
             if (nrow(tree1$cptable)==1){
               tree <- tree1}
             else {
               cventry <- which.min(tree1$cptable[, "xerror"])
               if (no.SE == 0){
                 cpcv <- tree1$cptable[cventry, "CP"]
                 tree <- prune(tree1, cp=cpcv)}
               else {
                 xerrorcv <- tree1$cptable[cventry, "xerror"]
                 sexerrorcv <- xerrorcv + tree1$cptable[cventry, "xstd"] * no.SE
                 cpcvse <- tree1$cptable[which.max(tree1$cptable[, "xerror"] <= sexerrorcv), "CP"]
                 tree <- prune(tree1, cp=cpcvse)}
         	}
	  }
         else {
             tree <- rpart(formula(paste(c("AdjustedTarget", Predictors),
                 collapse = "~")), data = newdata, subset = subs,
                 method = "anova", control = tree.control)
         }
	if(verbose) print(tree)

	## Estimate New Random Effects and Errors using LME
	# Get variables that identify the node for each observation
	newdata[,"nodeInd"] <- 0
	newdata[subs,"nodeInd"] <- tree$where
	# Fit linear model with nodes as predictors (we use the original target so likelihoods are comparable)
	# Check that the fitted tree has at least two nodes.
        if(min(tree$where)==max(tree$where)){
	    lmefit <- lme(formula(paste(c(toString(TargetName),1), collapse="~")), data=newdata, random=random, 
			subset=SubsetVector, method=method, control=lme.control, correlation=correlation)
	} else {
	    lmefit <- lme(formula(paste(c(toString(TargetName),"as.factor(nodeInd)"), collapse="~")), data=newdata, random=random, 
			subset=SubsetVector, method=method, control=lme.control, correlation=correlation)
        }
         adjtarg <- unique(cbind(tree$where, predict(lmefit, level=0)))
         tree$frame[adjtarg[,1],]$yval <- adjtarg[,2]

	if(verbose){
	    print(lmefit)
	    print(paste("Estimated Error Variance = ", lmefit$sigma))
	    print("Estimated Random Effects Variance = ")
	    print(as.matrix(lmefit$modelStruct$reStruct[[1]])*lmefit$sigma^2)
	}

	# Get the likelihood to check on convergence
	newlik <- logLik(lmefit)
	if(verbose) print(paste("Log likelihood: ", newlik))

	ContinueCondition <- (newlik-oldlik>ErrorTolerance & iterations < MaxIterations)
	oldlik <- newlik

	# Extract random effects to make the new adjusted target
	AllEffects <- lmefit$residuals[,1]-lmefit$residuals[,dim(lmefit$residuals)[2]]
	AdjustedTarget[subs] <- Target[subs] - AllEffects
    }

    residuals <- rep(NA, length=length(Target))
    residuals[subs] <- Target[subs]-predict(lmefit)
    attr(residuals, "label") <- NULL


   adjtarg <- unique(cbind(tree$where, predict(lmefit, level=0)))
   tree$frame[adjtarg[,1],]$yval <- adjtarg[,2]



    result <- list(Tree=tree, EffectModel=lmefit, RandomEffects=ranef(lmefit),
		BetweenMatrix=as.matrix(lmefit$modelStruct$reStruct[[1]])*lmefit$sigma^2,
		ErrorVariance=lmefit$sigma^2, data=data, logLik=newlik,
		IterationsUsed=iterations, Formula=formula, Random=random, Subset=subs,
		ErrorTolerance=ErrorTolerance, correlation=correlation,
		residuals=residuals, method=method, cv=cv, lme.control=lme.control, tree.control=tree.control)
    class(result) <- "REEMtree"

    return(result)
}



### Function: plot.REEMtree
# This function plots the tree part of a RE-EM tree,
# along with its text.
# Inputs:
#	x - a RE-EM tree object
#	text [TRUE] - should the text of the tree be included?
plot.REEMtree <- function(x, text=TRUE,...){
    plot(x$Tree)
    text(x$Tree)
}

### Function: print.REEMtree
# This function prints the RE-EM tree and the 
# associated covariances.
print.REEMtree <- function(x,...){
    print("*** RE-EM Tree ***")
    print(x$Tree)
    print("Estimated covariance matrix of random effects:")
    print(x$BetweenMatrix)
    print(paste("Estimated variance of errors:",x$ErrorVariance))
    print(paste("Log likelihood: ", x$logLik))
}

### Function: loglik.REEMtree
# This function computes the log likelihood of an estimated
# RE-EM tree.
logLik.REEMtree <- function(object,...){
    return(object$logLik)
}

### Function: predict.REEMtree
# This function predicts new observations for a RE-EM
# tree, given a dataset with columns of the same names.
# If the options EstimateRandomEffects is true, then
# this function will check the dataset for non-missing 
# values of the target with which to estimate the 
# random effects for observations; otherwise, it will
# use only the tree estimates, which sets the random
# effects to 0.
# NOTE: THIS ASSUMES ONLY RANDOM INTERCEPTS IN THE CASE
# WHERE RANDOM EFFECTS ARE ESTIMATED!
# Inputs:
# object - RE-EM tree
# newdata - new data set (with the same variable names)
# id [NULL] - a vector containing the identifiers for the groups 
# 	(needed only if random effects are estimated and newdata is not the data used for estimation)
# EstimateRandomEffects [TRUE] - should random effects be estimated from the given data?
# ... - other arguments to be passed through to RPART
# Output: a vector of predicted values, with one
#	for each element of the dataset
predict.REEMtree <- function(object, newdata, id=NULL, EstimateRandomEffects=TRUE,...){
    # Base predictions from the tree part
    treePrediction <- predict(object$Tree,newdata,...)

    # If we aren't estimating random effects, we 
    # just use the tree for prediction.
    if(!EstimateRandomEffects){
	return(treePrediction)
    }

    # Get the group identifiers if necessary
    if(is.null(id)){
	id <- object$EffectModel$groups
    }

    # Error-checking: the number of observations
    # in the dataset must match the sum of NumObs
    if(length(id) != dim(newdata)[1]){
	stop("number of observations in newdata does not match the length of the group identifiers")
    }

    ### Use the formula to get the target name
    TargetName <- object$Formula[[2]]
    # Remove the name of the data frame if necessary
    if(length(TargetName)>1) 
	TargetName <-TargetName[3]
    ActualTarget <- newdata[,toString(TargetName)]

    completePrediction <- treePrediction


    # Get the identities of the groups in the data
    # This will be slow - does LME have a faster way?
    uniqueID <- unique(id)

    # Get the random effects from the estimated RE-EM tree, in case there is overlap
    estRE <- ranef(object)
    
    for(i in 1:length(uniqueID)){
        # Identify the new group in the data
	nextID <- uniqueID[i]
	thisGroup <- id==nextID

	# If this group was in the original estimation, apply its random effect
	estEffect <- estRE[toString(uniqueID[i]),1]

	if(is.na(estEffect)){
	    # Check for non-missing target
	    nonMissing <- !is.na(ActualTarget[thisGroup])
	    numAvailable <- sum(nonMissing)

	    # If all the targets are missing, accept the
	    # tree prediction; otherwise, estimate
	    if(numAvailable>0) {
	    	R <- object$ErrorVariance * diag(numAvailable)
	    	D <- object$BetweenMatrix
	    	Z <- matrix(data=1,ncol=1, nrow=numAvailable)
    	    	W <- solve(R + Z %*% D %*% t(Z))
	    	effect <- D %*% t(Z) %*% W %*% 
			subset(ActualTarget[thisGroup] - treePrediction[thisGroup],subset=nonMissing)
	    	completePrediction[thisGroup] <- treePrediction[thisGroup]+effect
	    }
	} else {
	    completePrediction[thisGroup] <- treePrediction[thisGroup]+estEffect
	}

    }

    return(completePrediction)
}

### Function: ranef.REEMtree
# This function extracts the vector of estimated random effects
# from a RE-EM tree.  
ranef.REEMtree <- function(object,...){
    return(object$RandomEffects)
}

### Function: tree
# This function extracts the estimated tree from a RE-EM tree.
tree.REEMtree <- function(object,...){
    return(object$Tree)
}

### Function: tree
# This is the generic (S3-type) of function that returns the tree associated
# with an object (such as a RE-EM tree or a FE-EM tree).
tree <- function(object,...){
    UseMethod("tree")
}


### Function: is.REEMtree
# This function returns TRUE if an object has class REEMtree.
is.REEMtree <- function(object){
    return(class(object)=="REEMtree")
}

### Function: AutoCorrelationLRtest
# This function takes in a RE-EM tree and computes a likelihood ratio test 
# to check if there is autocorrelation (of the AR(1) type) in the errors).
# Note that this test keeps the tree structure fixed and re-estimates the 
# LME part only.  
# Inputs:
#	object - a RE-EM tree
#	newdata [NULL] - dataset on which the test is to be performed (the original dataset by default)
#	correlation [corAR1()] - type of autocorrelation to test
# Output: a list containing:
#     	loglik0 - the likelihood if 0 autocorrelation is assumed
#	loglikAR - the likelihood if autocorrelation is allowed
#	pvalue - the p-value from the chi-squared(1) distribution corresponding to logratio
AutoCorrelationLRtest <- function(object, newdata=NULL, correlation=corAR1()){
    # Extract the inputs we need from the RE-EM tree
    tree <- object$Tree
    formula <- object$Formula
    TargetName <- formula[[2]]
    # Remove the name of the data frame if necessary
    if(length(TargetName)>1) 
	TargetName <-TargetName[3]
    random <- object$Random
    method <- object$method
    lme.control <- object$lme.control

    if(is.null(newdata)) newdata <- object$data

    newdata[,"nodeInd"] <- tree$where

    # Fit the two linear models with nodes as predictors, with and without AR(1) correlation
    # Check that the fitted tree has at least two nodes.
    if(min(newdata$nodeInd)==max(newdata$nodeInd)){
	lme0 <- lme(formula(paste(c(toString(TargetName),1), collapse="~")), data=newdata, random=random, 
			method=method, control=lme.control, correlation=NULL)
	lmeAR <- lme(formula(paste(c(toString(TargetName),1), collapse="~")), data=newdata, random=random, 
			method=method, control=lme.control, correlation=correlation)
    } else {
	lme0 <- lme(formula(paste(c(toString(TargetName),"as.factor(nodeInd)"), collapse="~")), data=newdata, random=random, 
			method=method, control=lme.control, correlation=NULL)
	lmeAR <- lme(formula(paste(c(toString(TargetName),"as.factor(nodeInd)"), collapse="~")), data=newdata, random=random, 
			method=method, control=lme.control, correlation=correlation)
    }

    loglik0 <- logLik(lme0)
    loglikAR <- logLik(lmeAR)


    # Print out results
    print("*** Likelihood ratio test for autocorrelation ***")
    print(paste("Log likelihood for the hypothesis of no autocorrelation:", loglik0))
    print(paste("Log likelihood for the hypothesis of AR(1) autocorrelation:", loglikAR))
    print(paste("Test statistic:",-2*(loglik0-loglikAR)))
    print(paste("Asymptotic P-value:", pchisq(-2*(loglik0-loglikAR), 1, lower.tail=FALSE)))

    return(invisible(list(correlation=correlation, loglik0=loglik0,loglikAR=loglikAR, 
		pvalue=pchisq(-2*(loglik0-loglikAR), 1, lower.tail=FALSE))))
}

### Function: fitted
# This function replicated the fitted() function of LME (so that you don't have to call
# fitted(reem.result$EffectModel) instead).
fitted.REEMtree <- function(object, level=-1, asList=FALSE,...){
    if(level==-1){
	level <- dim(object$EffectModel$fitted)[2]-1
	print(level)
    } 
    return(fitted(object$EffectModel, level, asList))
}

### Function: residuals
# This function replicated the residuals() function of LME (so that you don't have to call
# residuals(reem.result$EffectModel) instead).
residuals.REEMtree <- function(object,level=-1, type="pearson", asList=FALSE,...){
    if(level==-1){
	level <- dim(object$EffectModel$fitted)[2]-1
	print(level)
    } 
    return(residuals(object$EffectModel, level, type,asList))
}








