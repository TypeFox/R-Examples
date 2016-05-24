`predict.oblique.tree` <- function(
	object,						#object of class "tree"
	newdata,					#data frame containing the values at which predictions are required. The 
	type = c("vector", "tree", "class", "where"),	#character string denoting whether the predictions are returned as a vector (default) or as a tree object. 
	eps = 1e-3,					#a lower bound for the probabilities, used if events of predicted probability zero occur in newdata when predicting a tree
	update.tree.predictions = FALSE,		#yval and yprob does not change with newdata though this can be changed
	...)
#DESCRIPTION
#	given object of class "tree", predictions of tree on training data or upon new test data is returned in form as requested		
#OUTPUT		type = "vector"		vector of predicted responses or, if the response is a factor, matrix of predicted class probabilities. This new object is obtained by dropping newdata down object. For factor predictors, if an observation contains a level not used to grow the tree, it is left at the deepest possible node and frame$yval or frame$yprob at that node is the prediction. 
#			type = "tree"		an object of class "tree" is returned with new values for frame$n and frame$dev. If newdata does not contain a column for the response in the formula the value of frame$dev will be NA, and if some values in the response are missing, the some of the deviances will be NA. 
#			type = "class"		for a classification tree, a factor of the predicted classes (that with highest posterior probability, with ties split randomly). 
#			type = "where"		the nodes the cases reach.
{
	#######################################################################
	impurity.predict <- function(
		yprob,
		newdata.actual.classes)
	{
		if (length(Y) == 0) {
			impurity <- 0
		} else {
			yprob[yprob == 0] <- max(0,eps)
			impurity <- -2*sum(log(yprob)*table(newdata.actual.classes))
		}

		return(impurity)
	}
	#######################################################################

	#get data
	type <- match.arg(type)

	#predict training data or new data?
	if (!missing(newdata)) {
		#######################################################################
		predict.recurse <- function(	X,		#data frame containing subsetted original dataset with model.frame() applied to it
						Y,
						node.name)	#running counter of the natural naming of nodes in a sequential manner
		#DESCRIPTION
		#	given data and the method with which node impurities are calculated, new examples are predicted to $n, $dev, $where and $y for the new data
		#OUTPUT		NOTHING
		{
#print(node.name)
#browser()

			#recursively predict tree
			object$frame$n[results.frame.row.counter] <<- length(Y)

			#do tree predictions change with new data?
			if (update.tree.predictions) {
				#tree structure changes with new data, update $dev, $yval and $yprob as well as $n
				class.probabilities <- table(Y) / length(Y)
				object$frame$dev[results.frame.row.counter] <<- node.impurity(	class.probabilities = class.probabilities,
												impurity.measure = split.impurity
										) * length(Y)
				object$frame$yval[results.frame.row.counter] <<- ylevels[which.is.max(class.probabilities)]
				object$frame$yprob[results.frame.row.counter,] <<- class.probabilities
			} else {
				#tree structure doesnt change, only update $dev as well as $n
				object$frame$dev[results.frame.row.counter] <<- impurity.predict(	yprob = object$frame$yprob[results.frame.row.counter,],
													newdata.actual.classes = Y)
			}


			#what is the current node?
			split.attribute <- object$frame$var[results.frame.row.counter]
			if (split.attribute == "<leaf>") {
				#its a leaf, stop recursing here
				observations.at.leaf <- names(results.where)%in%dimnames(X)[[1]]
				results.where[observations.at.leaf] <<- results.frame.row.counter
				results.y[observations.at.leaf] <<- object$frame$yval[results.frame.row.counter]	
				results.frame.row.counter <<- results.frame.row.counter + 1
			} else {
				#its an internal node, what type of split is it?
				if (split.attribute == "") {
					#its an oblique split
					coefs <- object$details[[results.frame.row.counter]]$coefficients
					coefs <- coefs[c(TRUE,coefs[-1]!=0)]
#browser()
					#treat if correctly, axis-parallel splits hidden as oblique splits are treated incorrectly
					if (length(coefs)==2) {
						#axis-parallel splits masquerading as oblique splits
						child.left.T.F <- X[,dimnames(X)[[2]] %in% names(coefs[2])] < -coefs[1]/coefs[2]
					} else {
						#oblique split with more than 1 continuous attribute
						child.left.T.F <- coefs[1] + X[,dimnames(X)[[2]] %in% names(coefs[-1])] %*% coefs[-1] < 0
					}
#print(child.left.T.F)
#browser()
				} else {
					#its not an oblique split, its univariate of some sort
					attribute.data <- X[,dimnames(X)[[2]] == split.attribute]
					attribute.type <- class(attribute.data)[1]	#ordered attributes have class "ordered" "factor" whilst factor attributes have class "factor"

					#apply the split
					if (attribute.type == "factor") {
						#categorical attribute
#browser()
						child.left.T.F <- 
							eval(	substitute(	data$foo%in%object$details[[results.frame.row.counter]]$variables,
										list(	foo=as.character(object$frame$var[[results.frame.row.counter]])
										)
								)
							)
					} else if (attribute.type == "ordered") {
						#ordered categorical attribute
#browser()
						child.left.T.F <- 
							eval(	substitute(	data$foo%in%object$details[[results.frame.row.counter]]$variables,
										list(	foo=as.character(object$frame$var[[results.frame.row.counter]])
										)
								)
							)
					} else if (attribute.type == "numeric") {
						#axis-parallel split on continuous attribute
						child.left.T.F <- attribute.data < object$details[[results.frame.row.counter]]
					}
				}

				#increment row counter
				results.frame.row.counter <<- results.frame.row.counter + 1

				#recurse on left child
				predict.recurse(	X = X[child.left.T.F,,drop=FALSE],
							Y = Y[child.left.T.F],
							node.name = 2 * node.name)

				#recurse on right child
				predict.recurse(	X = X[!child.left.T.F,,drop=FALSE],
							Y = Y[!child.left.T.F],
							node.name = 2 * node.name + 1)
			}
			#returns nothing, it changes information in different environment on the fly
		}
		#######################################################################

		#there is 'newdata', alter the tree object provided
		m <- model.frame(	formula = as.formula(eval(object$call$formula)),	#coerce newdata into useful form
					data = newdata)

		#extract the method 'tree' was built with
		split.impurity <- attributes(object)$split.impurity

		#extract responses and attributes
		Y <- model.extract(m, "response")				#Y 		<	>	contains the targets of each observation
		X <- tree:::tree.matrix(m)					#X		<	>	matrix of the predictors
		ylevels <- levels(Y)						#ylevels 	<	>	contains the levels of the targets

		#create data structures to hold data, predict with new data
		results.frame.row.counter <- 1					#frame row counter
		results.where <- vector("numeric",length(Y))
		names(results.where) <- dimnames(X)[[1]]
		results.y <- results.where

		#build tree recursively
		predict.recurse(	X=X,
					Y=Y,
					node.name=1)				#start at root node

		#update tree with new information
#		object$frame$n <- results.frame.n				#!!!!!!!!!!!!!!!!!!!!	n
#		object$frame$dev <- results.frame.dev				#!!!!!!!!!!!!!!!!!!!!	dev
		results.y.to.factor <- factor(results.y)
		levels(results.y.to.factor) <- levels(object$y)
		object$y <- results.y.to.factor
		object$where <- results.where
	}

	#return correct form of the result
	if (type=="vector") {
		#return vector of probabilities, construct predicted vector
		columns <- length(attr(object,"ylevel"))
		result.matrix <-	matrix(	0,
						nrow=length(object$where),
						ncol=columns,
						dimnames=	list(	names(object$where),
									attr(object,"ylevel")
								)
					)
		leaf.numbers <- unique(object$where)		#node row numbers as opposed to node names
		for (leaf.numbers.index in 1:length(leaf.numbers)) {
			result.matrix.subset <- object$where==leaf.numbers[leaf.numbers.index]
			result.matrix[result.matrix.subset,] <- 	matrix(	rep(	object$frame$yprob[leaf.numbers[leaf.numbers.index],],
											times=sum(result.matrix.subset)
										),
										ncol=columns,
										byrow=TRUE
									)
		}
		return (result.matrix)
	} else if (type=="tree") {
		#return tree
		return(object)
	} else if (type=="class") {
		#return class
		names(object$y) <- NULL		#conform with brian's code
		return(object$y)
	} else if (type=="where") {
		#return where
		return(object$where)
	}
}

