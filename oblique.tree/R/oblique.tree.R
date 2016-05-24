`oblique.tree` <- function (
	formula,							#fine	#how to build tree, dependent ~ independent
	data,								#fine	#data frame containing all relevant variables in formula
	subset,								#dunno	#subset of data to use
	control = tree.control(nobs, ...), 				#fine
	method = "recursive.partition", 				#fine	#method can be "model.frame"
	split.impurity = c(	"deviance",
				"gini"),					#impurity measure to use when growing tree
	model = FALSE,								#a model frame containing dependent to independent variables for fitting the tree, used in cv.oblique.tree()
	oblique.splits = c(	"only",					#fine	#consider only oblique splits on continuous attributes
				"on",					#fine	#consider both oblique splits and axis-parallel splits on continuous attributes
				"off"),					#fine	#consider only axis-parallel splits on continuous attributes
	variable.selection = c(	"none",					#fine	#consider full oblique splits
				"model.selection.aic",			#????	#consider oblique splits from ordinary logistic regression models stepAIC'ed with AIC
				"model.selection.bic",				#consider oblique splits from ordinary logistic regression models stepAIC'ed with BIC
				"lasso.aic",					#consider oblique splits from penalized logistic regression models regularized with the lasso (L1) where lambda is chosen by AIC
				"lasso.bic"),					#consider oblique splits from penalized logistic regression models regularized with the lasso (L1) where lambda is chosen by AIC
	...)
{
	#######################################################################
	node.split.categorical <- function(
		categorical.attribute,						#categorical attribute
		Y,								#targets
		attribute.index)
	#DESCRIPTION
	#	when given dataframe that had model.frame() applied on, so 1st column is class and others are all categorical, the best such split is found
	#OUTPUT		impurity			numeric impurity value of the best categorical split, Inf means don't use categorical splits
	#			optimal.variables	name of best categorical attribute to split upon
	#			optimal.split.left		characters explicitly naming optimal split for examples to the left
	#			optimal.split.right	characters explicitly naming optimal split for examples to the right
	#			child.left		dataframe with row names as example names containing at least the class attribute of examples to be split left
	#			child.right		dataframe with row names as example names containing at least the class attribute of examples to be split right
	#			details			list containing details about the split
	#				$variables		the variables used in the split
	{
		#initialise a storage container for the best split found over this categorical attribute
		split <- list(	impurity=Inf,
				child.left.T.F=NULL,			#vector or TRUE's and FALSE's indicating allocations of observations to child nodes when applying this split
				cutleft=NULL,
				cutright=NULL,
				details=NULL)

		#count all possible splits for this continuous attribute
		residual.factors <- unique(categorical.attribute)
		number.of.residual.factors <- length(residual.factors)

		#find best allocation of factors to child nodes

		#'number.of.residual.factors' may be 1, i.e. we cannot split on it
		if (number.of.residual.factors > 1) {
			#evaluate every possible superclass split on this categorical attribute
			for (superclass.index in 1:(2^(number.of.residual.factors-1)-1)) {
				#create the composition of the i^th superclass
				superclass.left <- generate.ith.superclass(	R = number.of.residual.factors,
										superclass.index = superclass.index)

				#evaluate allocations to child nodes for this split
				child.left.T.F <- categorical.attribute %in% residual.factors[as.logical(superclass.left)]
				number.in.child.left <- sum(child.left.T.F)
				number.in.child.right <- length(child.left.T.F) - number.in.child.left

				#only consider splits with >= mincut examples in either child
				if (number.in.child.left >= control$mincut && number.in.child.right >= control$mincut) {
					#evaluate average impurity of child nodes
					impurity <- 	(	node.impurity(	class.probabilities = table(Y[child.left.T.F]) / number.in.child.left,
										impurity.measure = split.impurity
								) * number.in.child.left +
								node.impurity(	class.probabilities = table(Y[!child.left.T.F]) / number.in.child.right,
										impurity.measure = split.impurity
								) * number.in.child.right
							) / length(Y)

					#check if resulting split is better than that found already, if so save result
					if (impurity < split$impurity) {
						split$impurity <- impurity
						split$cutleft <- paste(c(":",paste(xlevels[[attribute.index]][residual.factors[as.logical(superclass.left)]],collapse=",")),collapse="")
						split$cutright <- paste(c(":",paste(xlevels[[attribute.index]][residual.factors[!as.logical(superclass.left)]],collapse=",")),collapse="")
						split$child.left.T.F <- child.left.T.F
						split$details <- NULL
					}
				}
			}
		} #else there is only one class, i.e. we cannot split on this categorical attribute

		#return split
		return(split)
	}
	#######################################################################
	node.split.ordered.categorical <- function(
		categorical.attribute,						#ordered categorical attribute
		Y,								#targets
		attribute.index)
	#DESCRIPTION
	#	when given dataframe that had model.frame() applied on, so 1st column is class and others are all ordered categorical, the best such split is found
	#OUTPUT		impurity			numeric impurity value of the best ordered categorical split, Inf means don't use categorical splits
	#			optimal.variables	name of best categorical attribute to split upon
	#			optimal.split.left		characters explicitly naming optimal split for examples to the left
	#			optimal.split.right	characters explicitly naming optimal split for examples to the right
	#			child.left		dataframe with row names as example names containing at least the class attribute of examples to be split left
	#			child.right		dataframe with row names as example names containing at least the class attribute of examples to be split right
	#			details			list containing details about the split
	#				$variables		the variables used in the split
	{
		#initialise a storage container for the best split found over this categorical attribute
		split <- list(	impurity=Inf,
				child.left.T.F=NULL,			#vector or TRUE's and FALSE's indicating allocations of observations to child nodes when applying this split
				cutleft=NULL,
				cutright=NULL,
				details=NULL)

		#count all possible splits for this continuous attribute
		residual.factors <- unique(categorical.attribute)
		number.of.residual.factors <- length(residual.factors)

		#find best allocation of factors to child nodes

		#'number.of.residual.factors' may be 1, i.e. we cannot split on it
		if (number.of.residual.factors > 1) {
			#evaluate every possible superclass split on this categorical attribute
			for (superclass.index in 1:(number.of.residual.factors-1)) {
				#create the i^th superclass composition putting it in 'superclass.left'
				superclass.left <- numeric(number.of.residual.factors)
				superclass.left[1:superclass.index] <- 1

				#evaluate allocations to child nodes for this split
				child.left.T.F <- categorical.attribute %in% residual.factors[as.logical(superclass.left)]
				number.in.child.left <- sum(child.left.T.F)
				number.in.child.right <- length(child.left.T.F) - number.in.child.left

				#only consider splits with >= mincut examples in either child
				if (number.in.child.left >= control$mincut && number.in.child.right >= control$mincut) {
					#evaluate average impurity of child nodes
					impurity <- 	(	node.impurity(	class.probabilities = table(Y[child.left.T.F]) / number.in.child.left,
										impurity.measure = split.impurity
								) * number.in.child.left +
								node.impurity(	class.probabilities = table(Y[!child.left.T.F]) / number.in.child.right,
										impurity.measure = split.impurity
								) * number.in.child.right
							) / length(Y)
					#check if resulting split is better than that found already, if so save result
					if (impurity < split$impurity) {
						split$impurity <- impurity
						split$cutleft <- paste(c(":",paste(xlevels[[attribute.index]][residual.factors[as.logical(superclass.left)]],collapse=",")),collapse="")
						split$cutright <- paste(c(":",paste(xlevels[[attribute.index]][residual.factors[!as.logical(superclass.left)]],collapse=",")),collapse="")
						split$child.left.T.F <- child.left.T.F
						split$details <- NULL
					}
				}
			}
		} #else there is only one class, i.e. we cannot split on this categorical attribute

		#return split
		return(split)
	}
	#######################################################################
	node.split.continuous.axis.parallel <- function(
		continuous.attribute,						#continuous attribute
		Y)								#targets
	#DESCRIPTION
	#	when given dataframe that had model.frame() applied on, so 1st column is class and others are all continuous, the best axis parallel split is found
	#OUTPUT		impurity			numeric impurity value of the best categorical split
	#			optimal.variables	name of best categorical attribute to split upon
	#			optimal.split.left		characters explicitly naming optimal split for examples to the left
	#			optimal.split.right	characters explicitly naming optimal split for examples to the right
	#			child.left		dataframe with row names as example names containing at least the class attribute of examples to be split left
	#			child.right		dataframe with row names as example names containing at least the class attribute of examples to be split right
	{
		#initialise a storage container for the best split found over this continuous attribute
		split <- list(	impurity=Inf,
				child.left.T.F=NULL,			#vector or TRUE's and FALSE's indicating allocations of observations to child nodes when applying this split
				cutleft=NULL,
				cutright=NULL,
				details=NULL)

		#look at all possible splits for this continuous attribute
		unique.continous.data.points <- sort(unique(continuous.attribute))
		number.of.unique.splits <- length(unique.continous.data.points)

#write this like brian did and it'll speed up ur code significantly

		#if there is at least one split we can try, then apply it and look at the result
		if (number.of.unique.splits > 1) {
			#evaluate each of the 'number.of.unique.splits-1' splits on this continuous attribute
			for (split.index in 1:(number.of.unique.splits-1)) {
				cutoff.point <- (unique.continous.data.points[split.index] + unique.continous.data.points[split.index+1])/2

				#evaluate allocations to child nodes for this split
				child.left.T.F <- continuous.attribute < cutoff.point
				number.in.child.left <- sum(child.left.T.F)
				number.in.child.right <- length(child.left.T.F)-number.in.child.left

				#only consider splits with >= mincut examples in either child
				if (number.in.child.left >= control$mincut && number.in.child.right >= control$mincut) {
					#evaluate average impurity of child nodes
					impurity <- 	(	node.impurity(	class.probabilities = table(Y[child.left.T.F]) / number.in.child.left,
										impurity.measure = split.impurity
								) * number.in.child.left +
								node.impurity(	class.probabilities = table(Y[!child.left.T.F]) / number.in.child.right,
										impurity.measure = split.impurity
								) * number.in.child.right
							) / length(Y)
					#check if resulting split is better than that found already, if so save result
					if (impurity < split$impurity) {
						split$impurity <- impurity
						split$cutleft <- paste("<",round(cutoff.point,6),sep="")
						split$cutright <- paste(">",round(cutoff.point,6),sep="")
						split$child.left.T.F <- child.left.T.F
						split$details <- cutoff.point
					}
				}
			}
		} #else there is only one unique cutpoint, i.e. we cannot split on this continuous attribute

		return(split)
	}
	#######################################################################
	node.split.continuous.oblique <- function(
		all.continuous.attributes,					#all continuous attributes
		Y)								#targets
	#DESCRIPTION
	#	when given dataframe that had model.frame() applied on, so 1st column is class and others are all continuous, the best oblique split with our criteria is found
	#OUTPUT	impurity			numeric impurity value of the best continuous oblique split, Inf means don't use continuous oblique splits
	#		optimal.variables	name of best categorical attribute to split upon
	#		optimal.split.left		characters explicitly naming optimal split for examples to the left
	#		optimal.split.right	characters explicitly naming optimal split for examples to the right
	#		child.left		dataframe with row names as example names containing at least the class attribute of examples to be split left
	#		child.right		dataframe with row names as example names containing at least the class attribute of examples to be split right
	#		details			list containing details about the split
	#			$hyperplane.coefficients	hyperplane coefficients given accurately for later use
	#			$variables		the variables used in the hyperplane
	#			$superclass		the combination of classes g into sc1 (vs sc2) that result with best logistic regression decision boundary
	{
		#initialise a storage container for the best oblique split found
		best.split <- list(	impurity = Inf)

		#count number of possible oblique splits
		unique.levels <- levels(unique(Y))
		number.of.residual.factors <- length(unique.levels)

		#find best allocation of factors to child nodes (upon which linear boundaries can be found)
		#'number.of.residual.factors' may be 1, i.e. we cannot split on it
		if (number.of.residual.factors > 1) {
			#evaluate the full oblique splits resulting from every possible superclass allocation
			design.matrix <- 	cbind(	1,				#create design matrix for glm.fit() to use
							all.continuous.attributes
						)

			for (superclass.index in 1:(2^(number.of.residual.factors-1)-1)) {
				#create the composition of the i^th superclass
				superclass.left <- generate.ith.superclass(	R = number.of.residual.factors,
										superclass.index = superclass.index)

				#find the best full oblique split by extracting the linear boundary found by fitting an appropriate linear classifier to a smaller classification problem
				#create target information including the name of observations as well, will be useful for pruning purposes later on
				superclass.targets <- as.numeric(Y%in%levels(Y)[as.logical(superclass.left)])
				names(superclass.targets) <- row.names(all.continuous.attributes)
				logistic.regression.fit <- glm.fit(	x=design.matrix,
									y=superclass.targets,
									family=binomial())
#print(logistic.regression.fit$coef)
				#fitting to linearly separable data results in NA's which are 0, make them zero
				logistic.regression.fit$coefficients[is.na(logistic.regression.fit$coefficients)] <- 0

#glmpath$b.corrector[last.row,] is actually quite good!
				#evaluate allocation to child nodes for this oblique split
				child.left.T.F <- logistic.regression.fit$linear.predictors < 0
				number.in.child.left <- sum(child.left.T.F)
				number.in.child.right <- length(child.left.T.F) - number.in.child.left

#abline(a=-logistic.regression.fit$coefficients[1]/logistic.regression.fit$coefficients[3],b=-logistic.regression.fit$coefficients[2]/logistic.regression.fit$coefficients[3])

				#only consider splits with >= mincut examples in either child
				if (number.in.child.left >= control$mincut && number.in.child.right >= control$mincut) {
					#evaluate average impurity of child nodes
					impurity <- 	(	node.impurity(	class.probabilities = table(Y[child.left.T.F]) / number.in.child.left,
										impurity.measure = split.impurity
								) * number.in.child.left +
								node.impurity(	class.probabilities = table(Y[!child.left.T.F]) / number.in.child.right,
										impurity.measure = split.impurity
								) * number.in.child.right
							) / length(Y)
					#check if resulting oblique split is better than that found already, if so save result
					if (impurity < best.split$impurity) {
						#overwrite best split found so far with current split
						best.split <- oblique.split.writer(	coefficients = logistic.regression.fit$coefficients,
											impurity = impurity,
											child.left.T.F = child.left.T.F,
											superclass.composition = unique.levels[as.logical(superclass.left)])
					}
				}
			} #end for loop over all superclass oblique splits

##browser()

			#having considered all full oblique splits, has a legal split even been found?
			if (best.split$impurity != Inf) {
				#yes, a legal split was found

				#do we want to optimise the best oblique split to produce concise oblique splits during tree growth?
				if (variable.selection!="none") {
					#yes, we want to optimise the full oblique splits

					#find concise oblique split (in whatever manner we wish and update split$*, this requires that
					#	logistic.regression.fit$coefficients	contains updated coefficients of new split
					#	child.left.T.F				the updated T/F vector denoting to which child observations fall w.r.t. this new oblique split

#browser()

					#how should concise oblique splits be found?
					if (variable.selection %in% c("lasso.aic","lasso.bic")) {
						#use L1 regularisation to remove attributes
						logistic.regression.fit <- glmpath(	x=all.continuous.attributes,					#explanatory variables
											y=as.numeric(Y %in% best.split$details$superclass.composition),	#target variable
											family=binomial(),						#logistic regression
											lambda2=0,							#no L2 penalisation wanted
											standardize=FALSE,						#dont want to standardize explanatory variables to have unit variance
											max.steps=1000)

#
#plot(									glmnet(		x=all.continuous.attributes,					#explanatory variables
#											y=as.numeric(Y %in% best.split$details$superclass.composition),			#target variable
#											family="binomial",						#logistic regression
#											alpha=1,							#no L2 penalisation wanted
#											standardize=FALSE)
#)
#but i cant get hold of dev!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! cant automatically choose splits within family of subsplits

#print(logistic.regression.fit$lambda[1])
#
#	compare glmpath() and glm()
#
#	FORMULA=type~.
#	SAMP=sample(1:214,214)
#	DATA=fgl[SAMP,]
#	DATA$type=DATA$type%in%levels(DATA$type)[runif(6)>0.5]
#	X=tree:::tree.matrix(DATA[,1:9])
#	Y=as.numeric(DATA$type)
#	rbind(
#		glmpath(	x=X,
#				y=Y,
#				standardize=FALSE,
#				family=binomial(),
#				lambda2=0#
#		)$b.corrector,
#		glm(		formula=FORMULA,
#				data=DATA,
#				family=binomial()
#		)$coef
#	)


						#of the many L1 regularised splits, use split with lowest aic/bic that still produces a legal allocation of observations to child nodes
						#find row in logistic.regression.fit$b.corrector corresponding to penalised logistic regression fit with lowest aic or bic
#best.lambda.row <- which.min(logistic.regression.fit$aic)

#probably should rename active.set.lambda as its not lambda

						if (variable.selection == "lasso.aic") {
							#use aic info criterion
							active.set.lambda <- logistic.regression.fit$aic
						} else if (variable.selection == "lasso.bic") {
							#use bic info criterion
							active.set.lambda <- logistic.regression.fit$bic
						}
#browser()
						#collect out the lambda corresponding to changes in the active set, we only care about those ones, not everything found by the predictor-corrector method
						names(active.set.lambda) <- 1:length(active.set.lambda)

						#arriving at model with lowest aic/bic immediately, we may need to add attributes until a legal split is produced
						repeat {
							#check to see if legal splits are produced
							active.set.lambda.which.min <- which.min(active.set.lambda)
							best.active.set.lambda <- as.numeric(names(active.set.lambda.which.min))
							child.left.T.F <- predict(	object=logistic.regression.fit,
											newx=all.continuous.attributes,
											s=best.active.set.lambda,	#row number, that is, step number of the predictor-corrector algorithm
											type="link",
											mode="step") < 0		#evaluate predictions of penalised logistic regression model at points where attributes are introduced

							#evaluate allocations to child nodes for the split used
							number.in.child.left <- sum(child.left.T.F)
							number.in.child.right <- length(child.left.T.F) - number.in.child.left


###this doesnt allow lasso to produce a stump, i.e. no split, it only accepts legal splits
							#check to see if split is legal
							if (number.in.child.left >= control$mincut && number.in.child.right >= control$mincut) {
								#legal split produced, use this split
								break
							} else {
								#legal split NOT produced, add an attribute
								active.set.lambda <- active.set.lambda[-active.set.lambda.which.min]
							}
						} #we have obtained the most concise oblique split that produces legal allocations of observations to child nodes

						#package glmpath object with $coefficients passing information along about the hyperplane
						coefs <- logistic.regression.fit$b.corrector[best.active.set.lambda,]
						logistic.regression.fit$coefficients <- coefs[c(TRUE,coefs[-1]!=0)]			#store coefficients to logistic.regression.fit just like it is done for glm()s
					} else if (variable.selection %in% c("model.selection.aic","model.selection.bic")) {
						#use aic/bic to remove attributes

						#write correct value to stepaic.k which controls whether we are using aic or bic
						if (variable.selection == "model.selection.aic") {
							#use aic info criterion
							stepaic.k <- 2
						} else if (variable.selection == "model.selection.bic") {
							#use bic info criterion
							stepaic.k <- log(length(Y))
						}

						#only apply subset-selection to the best split found, recreate classification data
						design.matrix <- data.frame(	targets = as.numeric(Y %in% best.split$details$superclass.composition),
										all.continuous.attributes)

						#refit best oblique split found
						best.split$logistic.regression.fit <- glm(	formula = targets~.,
												family = binomial(),
												data = design.matrix)

						#optimise logistic regression model one step at a time by application of step() finding the model with lowest AIC or the smallest model still producing legal splits
						repeat {
							#one-step lookahead, current model in best.split$logistic.regression.fit, new model in logistic.regression.fit
							logistic.regression.fit <- step(	object = best.split$logistic.regression.fit,
												trace = FALSE,
												steps = 1,			#one-step lookahead
												k = stepaic.k)
							#defaults to "both", do stepwise search, NOT backwards OR forwards
#print("repeating step()\n")

							#evaluate allocations to child nodes for this split
							child.left.T.F <- logistic.regression.fit$linear.predictors < 0

							#evaluate allocations to child nodes for the split used
							number.in.child.left <- sum(child.left.T.F)
							number.in.child.right <- length(child.left.T.F) - number.in.child.left

							#check to see if split has been updated
							if (length(best.split$logistic.regression.fit$coefficients) != length(logistic.regression.fit$coefficients)) {
								#split updated, does it produce legal splits?

								#only consider splits with >= mincut examples in either child
								if (number.in.child.left >= control$mincut && number.in.child.right >= control$mincut) {
									#legal split produced, continue search for next legal split
									best.split$logistic.regression.fit <- logistic.regression.fit
								} else {
									#legal split NOT produced, previous split was best
									logistic.regression.fit <- best.split$logistic.regression.fit

									#evaluate allocations to child nodes for this split
									child.left.T.F <- logistic.regression.fit$linear.predictors < 0

									#evaluate allocations to child nodes for the split used
									number.in.child.left <- sum(child.left.T.F)
									number.in.child.right <- length(child.left.T.F) - number.in.child.left
									break
								}
							} else {
								#split not updated, use split$logistic.regression.fit, which is just saved in logistic.regression.fit
								break
							}
						} #end of repeat, should have the simplest legal split based on this allocation of superclasses
					} #best full oblique split has been optimised to find best concise oblique split

					#evaluate average impurity of child nodes
					impurity <- 	(	node.impurity(	class.probabilities = table(Y[child.left.T.F]) / number.in.child.left,
										impurity.measure = split.impurity
								) * number.in.child.left +
								node.impurity(	class.probabilities = table(Y[!child.left.T.F]) / number.in.child.right,
										impurity.measure = split.impurity
								) * number.in.child.right
							) / length(Y)

					#we are definitely choosing this oblique split whatever its value of impurity
					best.split <- oblique.split.writer(	coefficients = logistic.regression.fit$coefficients,
										impurity = impurity,
										child.left.T.F = child.left.T.F,
										superclass.composition = best.split$details$superclass.composition)
				} #else we only consider full oblique splits, job done, no need to find concise oblique splits
			} #else no legal split was found, job done, no split to even try optimising

#			#make sure splits consider have >= mincut examples in both child nodes
#			if (sum(best.split$child.left.T.F) < control$mincut || length(Y)-number.in.child.left < control$mincut) {
#				#split not legal, return no oblique split
#				best.split$impurity <- Inf
#			}


		} #else there is only one factor, no oblique split can be considered, job done

		#passes best.split along as well
		return(best.split)
	}
	#######################################################################
	node.split.best <- function(
		X,					#all attributes
		Y)					#targets
	#DESCRIPTION
	#	when given dataframe that had model.frame() applied on, so 1st column is class and others are all continuous, the best possible split is found over all p attributes
	#OUTPUT		impurity			numeric impurity value of the split with least impurity measure
	#			optimal.variables	name of attribute resulting with least impurity
	#			optimal.split.left		characters explicitly naming optimal split for examples to the left
	#			optimal.split.right	characters explicitly naming optimal split for examples to the right
	#			child.left		dataframe with row names as example names containing at least the class attribute of examples to be split left
	#			child.right		dataframe with row names as example names containing at least the class attribute of examples to be split right
	#			details			list containing details about the split
	#				$hyperplane.coefficients	hyperplane coefficients given accurately for later use, only exists for oblique splits
	#				$variables		the variables used in the hyperplane or categorical/ordered categorical splits
	#				$superclass		the combination of classes g into sc1 (vs sc2) that result with best logistic regression decision boundary, only exists for oblique splits
	{
#WHAT IF CATEGORIC AND ORDERED CATEGORIC VARIABLES ARE ALREADY PERFECTLY SAME CLASS?
		#initialise a storage container for the best split found so far
		#each of the node.split.* functions finds best split of its class overwriting best.split if impurity is strictly better, they should contain the following info
		#best.split <- list(	impurity=Inf,
		#			child.left.T.F=NULL,			#vector or TRUE's and FALSE's indicating allocations of observations to child nodes when applying this split
		#			cutleft=NULL,
		#			cutright=NULL,
		#			details=NULL,
		#			variable=NULL)
		best.split <- list(	impurity = Inf)

		#find best univariate split on both continuous and categorical attributes
		for (attribute.index in 1:dim(X)[2]) {
			#initialise 'split'
			split <- list(	impurity = Inf)

			#find the best split upon this attribute
			if (attributes(Terms)$dataClasses[attribute.index+1] == "factor") {
				#find best categorical split on X[,attribute.index]
				split <- node.split.categorical(		categorical.attribute = X[,attribute.index],
										Y = Y,
										attribute.index)
			} else if (attributes(Terms)$dataClasses[attribute.index+1] == "ordered") {
				#find best ordered categorical split on X[,attribute.index]
				split <- node.split.ordered.categorical(	categorical.attribute = X[,attribute.index],
										Y = Y,
										attribute.index)
			} else if (oblique.splits!="only" && attributes(Terms)$dataClasses[attribute.index+1] == "numeric") {
				#find best continuous axis-parallel split on X[,attribute.index]
				split <- node.split.continuous.axis.parallel(	continuous.attribute = X[,attribute.index],
										Y = Y)
			}

			#check if resulting split (whatever it may be) is better than that found already, if so save result
			if (split$impurity < best.split$impurity) {
				best.split <- split
				best.split$variable <- attributes(X)$dimnames[[2]][attribute.index]
			}
		}

		if (oblique.splits != "off") {
			#find best linear-multivariate (oblique) split on continuous attributes
			continuous.attributes.T.F <- attributes(Terms)$dataClasses[-1]=="numeric"
			if (sum(continuous.attributes.T.F) > 1) {
				#find best oblique split upon all continuous variables
				split <- node.split.continuous.oblique(		all.continuous.attributes = X[,continuous.attributes.T.F],
										Y = Y)
#abline(a=-logistic.regression.fit$coefficients[1]/logistic.regression.fit$coefficients[3],b=-logistic.regression.fit$coefficients[2]/logistic.regression.fit$coefficients[3])

				#check if resulting oblique split is better than that found already, if so save result
				if (split$impurity < best.split$impurity) {
					best.split <- split
				}
			}
		} #otherwise, dont consider oblique splits
############# browser()

		#return best split criterion
		return(best.split)
	}
	#######################################################################
	node.recurse <- function(
		X,								#matrix containing data from attributes of the 'nobs' observations
		Y,								#vector of targets for each observation
		node.name)							#running counter of the natural naming of nodes in a sequential manner
	#DESCRIPTION
	#OUTPUT		alters information in parent frame on the fly
	{
		#fill in basic information
		results.frame.row.names[results.frame.row.counter] <<- node.name
		results.frame.n[results.frame.row.counter] <<- length(Y)
			class.probabilities <- table(Y) / length(Y)
			dev <- 	node.impurity(	class.probabilities = class.probabilities,
						impurity.measure = split.impurity
				) * length(Y)
#
#interesting, i thought it was always "deviance"
#
		results.frame.dev[results.frame.row.counter] <<- dev
			yval <- ylevels[which.is.max(class.probabilities)]
		results.frame.yval[results.frame.row.counter] <<- yval
		results.frame.yprob[results.frame.row.counter,] <<- class.probabilities

		#decide whether current node should be a leaf or not in order to fill in rest of the node info
		if (length(levels(factor(Y))) == 1 || length(Y) < control$minsize || control$mindev * results.frame.dev[1] > dev) {
			#perfectly classified || fewer than min.size examples in this node || sufficiently pure, return leaf
			results.frame.var[results.frame.row.counter] <<- "<leaf>"
			#cutleft and cutright defaulted to "" so nothing to change
			observations.of.interest <- names(results.where) %in% row.names(X)
			results.where[observations.of.interest] <<- results.frame.row.counter	#node.name
			results.y[observations.of.interest] <<- yval
			results.frame.row.counter <<- results.frame.row.counter + 1
		} else {
			#we have more than one class and a sufficiently impure node, find best binary split
			best.split <- node.split.best(	X = X,
							Y = Y)
			#has a legal split been found?
			if (best.split$impurity == Inf) {
				#no, make current node a leaf
				results.frame.var[results.frame.row.counter] <<- "<leaf>"
				#cutleft and cutright already "" so nothing to change
				observations.of.interest <- names(results.where)%in%row.names(X)
				results.where[observations.of.interest] <<- results.frame.row.counter	#node.name
				results.y[observations.of.interest] <<- yval
				results.frame.row.counter <<- results.frame.row.counter + 1
			} else {
				results.frame.var[results.frame.row.counter] <<- best.split$variable
				cutleft[results.frame.row.counter] <<- best.split$cutleft
				cutright[results.frame.row.counter] <<- best.split$cutright
				details[[results.frame.row.counter]] <<- best.split$details
				results.frame.row.counter <<- results.frame.row.counter + 1

				#recurse on left child
				node.recurse(	X = X[best.split$child.left.T.F,,drop=FALSE],
						Y = Y[best.split$child.left.T.F],
						node.name = 2 * node.name)

				#recurse on right child
				node.recurse(	X = X[!best.split$child.left.T.F,,drop=FALSE],
						Y = Y[!best.split$child.left.T.F],
						node.name = 2 * node.name + 1)
			}
		}
	}
	#######################################################################

	#make 'm' the model frame
	if (is.data.frame(model)) {
		#if 'model' present, use it
		m <- model
		model <- FALSE
	} else {
		#if 'model' not, create appropriate model frame
		m <- match.call(expand.dots = FALSE)
		m$control <- m$method <- m$split.impurity <- m$model <- m$oblique.splits <- m$variable.selection <- m$... <- NULL
		m[[1]] <- as.name("model.frame.default")
		m <- eval.parent(m)
		if (method == "model.frame")
			return(m)
	}

	#extract correct arguments
	split.impurity <- match.arg(split.impurity)				#impurity measure to compare splits
	oblique.splits <- match.arg(oblique.splits)				#impurity measure to compare splits
	variable.selection <- match.arg(variable.selection)			#whether full logistic regression models or aic/bic optimised models should be used

	#check if interaction terms present
	Terms <- attr(m, "terms")						#Terms		<	>	contains the terms of the model frame
	if (any(attr(Terms, "order") > 1))
  		stop("trees cannot handle interaction terms")

	#check if multiple responses present
	Y <- model.extract(m, "response")					#Y 		<	>	contains the targets of each observation
	if (is.matrix(Y) && ncol(Y) > 1)
		stop("trees cannot handle multiple responses")
	ylevels <- levels(Y)							#ylevels 	<	>	contains the levels of the targets

	#if there are unknown targets, stop the function
	if (any(is.na(Y)))
		stop("there are unknown responses")

	#ensure offsets not provided for classification trees
	offset <- attr(Terms, "offset")						#offset		<	>	contains any offset terms that might be used
	if (!is.null(offset)) {
		#no offset given
		stop("offset not implemented for classification trees")
#		offset <- m[[offset]]
#		Y <- Y - offset
	}

	#simplify model frame and associated data to optimize latter calculations
	X <- tree:::tree.matrix(m)						#X		<	>	matrix of the predictors
	xlevels <- attr(X, "column.levels")					#xlevels	<	>	the levels of factors in each attribute
	if (is.null(xlevels)) {
		xlevels <- rep(list(NULL), ncol(X))
		names(xlevels) <- dimnames(X)[[2]]
	}

	#check if there are actually any observations and if that number satisfies our controls
	nobs <- length(Y)							#nobs		<	>	number of observations (including NA's)
	if (nobs == 0)
		stop("no observations from which to fit a model")
	if (!is.null(control$nobs) && control$nobs < nobs) {
		stop("control$nobs < number of observations in data")
	}

	#make containers to store information from oblique tree that will be grown recursively
	results.frame.row.counter <- 1						#pointer to row of results.frame to print data
	results.frame.row.names <- vector("numeric",control$nmax)
	results.frame.var <- vector("character",control$nmax)
	results.frame.n <- vector("numeric",control$nmax)
	results.frame.dev <- vector("numeric",control$nmax)
	results.frame.yval <- vector("character",control$nmax)
	cutleft <- vector("character",control$nmax)
	cutright <- vector("character",control$nmax)
	results.frame.yprob <- matrix(0,nrow=control$nmax,ncol=length(ylevels),dimnames=list(NULL,ylevels))
	details <- vector("list",control$nmax)
	results.where <- vector("numeric",nobs)
	names(results.where) <- row.names(m)
	results.y <- results.where						#they have the same structure, recycle it
	results.weights <- rep(1,times=nobs)

	#grow a oblique tree in R
	options("warn"=-1)							#use glm() silently
	node.recurse(
		X =		X,
		Y =		Y,
		node.name =	1)						#start at root node
	options("warn"=0)

	#housekeeping
	frame <- data.frame(	var = 		factor(	results.frame.var[1:(results.frame.row.counter-1)],
							levels=c("<leaf>","",dimnames(X)[[2]])
						), #levels
				n = 		results.frame.n[1:(results.frame.row.counter-1)],
				dev = 		results.frame.dev[1:(results.frame.row.counter-1)],
				yval =		factor(	results.frame.yval[1:(results.frame.row.counter-1)],
					 	levels=ylevels),
				row.names = 	as.integer(results.frame.row.names[1:(results.frame.row.counter-1)])
		)	#given correct factors
	frame$splits <- cbind(cutleft,cutright)[1:(results.frame.row.counter-1),]
	frame$yprob <- results.frame.yprob[1:(results.frame.row.counter-1),]

	results <- 	list(	frame =		frame,
				details =	details[1:(results.frame.row.counter-1)]
			)
	results$where <- 	results.where					#which leaf examples fall into
	results$terms <- 	Terms
	results$call <- 	match.call()					#call
	results$y <- 		factor(results.y,levels=ylevels)		#predicted values
	results$weights <- 	results.weights					#weights for examples

	#coerce final result into object of class c("oblique.tree","tree") to make use of functions in the tree library
	class(results) <- c("oblique.tree","tree")
	attr(results, "xlevels") <- xlevels
	attr(results,"ylevels") <- ylevels
	attr(results,"split.impurity") <- split.impurity
	return(results)
}

