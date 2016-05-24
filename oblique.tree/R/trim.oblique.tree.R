`trim.oblique.tree` <- function(
	tree,						#object of class "tree"
#	k = NULL,					#vector containing cost complexities to prune to, if empty, summarizes entire tree sequence, if numeric vector given, tree sequence of trees pruned to such a cost complexities is summarized and if just a single value, returns tree pruned to such a cost complexity
	best = NULL,					#?
	newdata,					#new data with which to be used to update tree$frame
	trim.impurity = c("deviance","misclass"),	#impurity measure to use with which to compare subtrees for pruning
	trim.depth = c("partial","complete"),		#controls the type of splits considered
	eps = 1e-3)
#DESCRIPTION
#	given object of class "tree", impurity measure to with which compare subtrees, will return summary of tree sequence if cost complexity k is NOT specified or it a numeric vector (greater than length 1) is specified, if of length 1, will return optimal subtree pruned to that cost complexity
#OUTPUT		impurity			numeric impurity value of subtree with such a composition of examples
{
	#######################################################################
	trim.once <- function(
		tree,						#object of class "tree" to be pruned
		leaf.indicator)
###!!! dont actually need to pass this, can just figure it out again but it wastes a little calculation
	#DESCRIPTION
	#given a tree, the next optimal tree in the optimal tree sequence is found and returned along with all g(t) for all subtrees of the tree
	#OUTPUT		$tree=tree		object of class "tree" that resulting from pruning to the next optimal cost complexity value
	#			$g.of.t			the value of g(t) for all subtrees of the tree
	{
#print("choose nodes to prune or trim\n");
#browser()
#original.tree<-tree
		#get row numbers of internal nodes
		internal.nodes.row.numbers <- which(!leaf.indicator)

		#and the names of the leaves
		leaf.node.names <- as.numeric(row.names(tree$frame)[leaf.indicator])

		#create containers to hold information about subtrees (their leaves and g(t))
		subtree.node.names <- subtree.leaf.names <- trimmed.trees <- vector("list",length(internal.nodes.row.numbers))
		g.of.t.for.trims <- g.of.t.for.prunes <- rep(Inf,length(internal.nodes.row.numbers))	#everything defaulted to Inf, i.e. no trim is possible

		#compute g.of.t
		for (subtree.index in 1:length(internal.nodes.row.numbers)) {
			#1 of 4 things can be true
			#		|	oblique			not oblique
			#---------------+------------------------------------------
			#complete	|	R(T_t), R(T_trim)	R(T_t), R(t)
			#partial	|	R(T_t), R(T_trim)	nothing
			subtree.row.number <- internal.nodes.row.numbers[subtree.index]
#	print(subtree.row.number)
#	browser()
			if (trim.depth == "complete" || length(tree$details[[subtree.row.number]]) == 3) {
				#either "complete" (all nodes considered for trimming/pruning) or current node applies an oblique split, compute R.of.T.t either way
				#		|	oblique			not oblique
				#---------------+------------------------------------------
				#complete	|	R(T_t), R(T_trim)	R(T_t), R(t)
				#partial	|	R(T_t), R(T_trim)

				#internal node applies an oblique split
				subtree.root.name <- as.numeric(row.names(tree$frame)[subtree.row.number])

				#find last node for subtree T_'subtree.root.name'
				last.node.row.number <-	last.node.of.subtree(	tree = tree,
										subtree.root.name = subtree.root.name
							)

				#last.node.in.subtree.row.number now contains the number of the last node in subtree of interest
				subtree.node.names[[subtree.index]] <- as.numeric(row.names(tree$frame)[subtree.row.number:last.node.row.number])
				subtree.leaf.names[[subtree.index]] <- subtree.node.names[[subtree.index]][-1][subtree.node.names[[subtree.index]][-1] %in% leaf.node.names]

				#create a vector of T/F to indicate rows of tree$frame that are leaves from this subtree
				subtree.leaves.rows.T.F <- row.names(tree$frame) %in% subtree.leaf.names[[subtree.index]]
#browser()
				#evaluate measure of misfit for subtree T.t
				R.of.T.t <- tree.impurity(	yprob = tree$frame$yprob[subtree.leaves.rows.T.F,,drop=FALSE],
								number.of.observations.at.leaves = tree$frame$n[subtree.leaves.rows.T.F],
								leaf.classes = tree$frame$yval[subtree.leaves.rows.T.F],
								impurity.measure = trim.impurity)

				#evaluate complexity penalty for subtree T.t
				penalty.of.R.of.T.t <-	oblique.tree.complexity(	tree = tree,
											subtree.internal.node.names = subtree.node.names[[subtree.index]][!(subtree.node.names[[subtree.index]] %in% leaf.node.names)]
							)

				#what we trim to depends on what type of split is applied at the current internal node
				if (length(tree$details[[subtree.row.number]]) == 3) {
#	print(c("oblique",subtree.index))
#	browser()

					#this node applies an oblique split, find R.of.T.t.trimmed
					#		|	oblique			not oblique
					#---------------+------------------------------------------
					#complete	|	R(T_t), R(T_trim)
					#partial	|	R(T_t), R(T_trim)

					#evaluate measure of misfit for trimmed subtree T.t.trimmed
					trimmed.trees[[subtree.index]] <- tree
					trimmed.trees[[subtree.index]]$details[[subtree.row.number]]$coefficients <- 	tree$details[[subtree.row.number]]$trim.sequence[	dim(tree$details[[subtree.row.number]]$trim.sequence)[1],
															]

					#look at what happens to observations at node of interest when node of interest is trimmed
					subtree.leaf.row.names <- which(row.names(tree$frame) %in% subtree.leaf.names[[subtree.index]])
					observations.in.subtree.of.interest <- tree$where %in% subtree.leaf.row.names

					trimmed.trees[[subtree.index]] <- predict.oblique.tree(	object = trimmed.trees[[subtree.index]],					#trimmed tree
												newdata = original.data[observations.in.subtree.of.interest,],	#only focusing on subset of data that falls into node of interest
												type = "tree",
												eps = eps,
												update.tree.predictions = TRUE)							#also updating $yprob and $yval rather than keeping them static and only changing $n and $dev, i.e. we preserve the tree structure trimming only one oblique split and we observe the changes in the final tree

#	#
#	#compare trimmed.trees[] and tree[]
#	#
#	#compare frames where there are observations
#	cbind(trimmed.trees[[subtree.index]]$frame[trimmed.trees[[subtree.index]]$frame$n != 0,-5],tree$frame[trimmed.trees[[subtree.index]]$frame$n != 0,-5])
#	as.matrix(trimmed.trees[[subtree.index]]$frame[trimmed.trees[[subtree.index]]$frame$n != 0,-5]) == as.matrix(tree$frame[trimmed.trees[[subtree.index]]$frame$n != 0,-5])
#
#	#compare entire frames
#	cbind(trimmed.trees[[subtree.index]]$frame[-5],tree$frame[-5])
#	trimmed.trees[[subtree.index]]$frame[,-5]-tree$frame[-5]
#	#nodes without observations in
#	as.numeric(row.names(tree$frame)[trimmed.trees[[subtree.index]]$frame$n != 0])

					#having "updated" tree growth by trimming an oblique split, NaN's may appear when there are 0 observations at a node, NaN's mess up tree.impurity() and $dev calculations
					#avoid by simply overwriting NaNs with 0
					presence.of.NaNs <- trimmed.trees[[subtree.index]]$frame$n == 0
					trimmed.trees[[subtree.index]]$frame$dev[presence.of.NaNs] <- 0
					trimmed.trees[[subtree.index]]$frame$yprob[presence.of.NaNs,] <- 0
					trimmed.trees[[subtree.index]]$frame$yval[presence.of.NaNs] <- tree$frame$yval[presence.of.NaNs]		#arbitrarily predict everything to class 1,

					R.of.T.t.trimmed <- tree.impurity(	yprob = trimmed.trees[[subtree.index]]$frame$yprob[subtree.leaves.rows.T.F,,drop=FALSE],
										number.of.observations.at.leaves = trimmed.trees[[subtree.index]]$frame$n[subtree.leaves.rows.T.F],
										leaf.classes = trimmed.trees[[subtree.index]]$frame$yval[subtree.leaves.rows.T.F],
										impurity.measure = trim.impurity)

					#evaluate complexity penalty for subtree T.t.trimmed
					penalty.of.R.of.T.t.trimmed <-	oblique.tree.complexity(	tree = trimmed.trees[[subtree.index]],
													subtree.internal.node.names = subtree.node.names[[subtree.index]][!(subtree.node.names[[subtree.index]] %in% leaf.node.names)]
									)

#	print(	c(	"R.of.T.t.trimmed",R.of.T.t.trimmed,
#			"R.of.T.t",R.of.T.t,
#			"penalty.of.R.of.T.t",penalty.of.R.of.T.t,
#			"penalty.of.R.of.T.t.trimmed",penalty.of.R.of.T.t.trimmed
#		)
#	)
#	browser()
					#evaluate g.of.t.for.trims
					g.of.t.for.trims[subtree.index] <- (R.of.T.t.trimmed - R.of.T.t) / (penalty.of.R.of.T.t - penalty.of.R.of.T.t.trimmed)

#print(c(penalty.of.R.of.T.t,penalty.of.R.of.T.t.trimmed))
#browser()
				} else {
#print(c("axis-parallel",subtree.index))
					#trim.depth is "complete" and current node is not oblique, axis-parallel split at this internal node "trimmed" to a leaf
					#		|	oblique			not oblique
					#---------------+------------------------------------------
					#complete	|				R(T_t), R(t)
					#partial	|

					#evaluate measure of misfit for pruned node t

					R.of.t <- tree.impurity(	yprob = tree$frame$yprob[subtree.row.number,,drop=FALSE],
									number.of.observations.at.leaves = tree$frame$n[subtree.row.number],
									leaf.classes = tree$frame$yval[subtree.row.number],
									impurity.measure = trim.impurity)

#	browser()
#	print(	c(	"R.of.t",R.of.t,
#			"R.of.T.t",R.of.T.t,
#			"penalty.of.R.of.T.t",penalty.of.R.of.T.t
#		)
#	)

					#evaluate g.of.t.for.prunes
					g.of.t.for.prunes[subtree.index] <- (R.of.t - R.of.T.t) / (penalty.of.R.of.T.t - 1)

				}
#g.of.t.for.trims <- g.of.t.for.prunes


			} # ELSE trim.depth is "partial" and current node is not oblique, skip to next internal node
			#		|	oblique			not oblique
			#---------------+------------------------------------------
			#complete	|
			#partial	|				nothing
		}


		#nodes with min(g.of.t) need to be trimmed
		min.g.of.t <- min(c(g.of.t.for.trims,g.of.t.for.prunes))

		#pruning trumps trimming, it pruning is does just as well, then prune and dont trim
		if (min.g.of.t %in% g.of.t.for.prunes) {
			#pruning is best or does just as well

			#prune nodes that need to be pruned
			nodes.to.prune.T.F <- g.of.t.for.prunes == min.g.of.t
			if (sum(nodes.to.prune.T.F) > 0) {
#print(g.of.t.for.prunes)
#print("look at tree before pruning\n");
#x11()
#par(mfrow=c(1,2))
#	subtree.leaves <- tree$frame$var == "<leaf>"
#	print(	oblique.tree:::tree.impurity(
#			yprob                                   = tree$frame$yprob[subtree.leaves,,drop=FALSE],
#			number.of.observations.at.leaves        = tree$frame$n[subtree.leaves],
#			leaf.classes                            = tree$frame$yval[subtree.leaves],
#			impurity.measure                        = "misclass")
#	)
#	print(sum(tree$y!=Pima.tr$type))
#before<-tree
#plot(before,type="u");try(text(before))
#browser()
				#node(s) needs pruning, alter	$frame and row.names($frame) to reflect the updated tree,
				#				$where to the new locations using the old node numbering system
				#				$y with the correct new predictions
				#whilst keep a running idea of which nodes have been pruned to allow removal of unused entries in mat and $details in one go
#print(c("prune",internal.nodes.row.numbers[nodes.to.prune.T.F]))
				tree <- oblique.tree.prune.nodes(	tree = tree,
									list.of.node.names.to.prune = subtree.node.names[nodes.to.prune.T.F])
			}
		} else {
#print(g.of.t.for.trims)
#print(g.of.t.for.prunes)
#print(paste("min.g.of.t = ",min.g.of.t))
#print("look at tree before trimming\n");
#x11()
#par(mfrow=c(1,2))
#	subtree.leaves <- tree$frame$var == "<leaf>"
#	print(	oblique.tree:::tree.impurity(
#			yprob                                   = tree$frame$yprob[subtree.leaves,,drop=FALSE],
#			number.of.observations.at.leaves        = tree$frame$n[subtree.leaves],
#			leaf.classes                            = tree$frame$yval[subtree.leaves],
#			impurity.measure                        = "misclass")
#	)
#	print(sum(tree$y!=Pima.tr$type))
#before<-tree
#plot(before,type="u");try(text(before))
#browser()

			#pruning does not do better, so we trim
			#but we only ever trim one node at a time to simplify matters, trimming too many changes tree predictions
			nodes.to.trim.T.F <- g.of.t.for.trims == min.g.of.t
			number.of.nodes.to.trim <- sum(nodes.to.trim.T.F)
#		if (number.of.nodes.to.trim > 0) {
			#node(s) needs to be trimmed, alter	$frame	reflecting the updated tree			##and row.names($frame)
			#					$where	new locations 					##using the old node numbering system
			#					$y	new predictions
#			list.of.node.names.to.trim = subtree.node.names[nodes.to.trim.T.F]
#			list.of.trimmed.trees = trimmed.trees[nodes.to.trim.T.F]
			list.of.node.names.to.trim = subtree.node.names[nodes.to.trim.T.F][number.of.nodes.to.trim]	#trim one node at a time, trim latter nodes first as earlier nodes will still have data to look at them again while latter nodes might not once earlier nodes are trimmed
			list.of.trimmed.trees = trimmed.trees[nodes.to.trim.T.F][number.of.nodes.to.trim]
#print(c("trim",internal.nodes.row.numbers[nodes.to.trim.T.F]))

#make sure tree trimming get axis-parallel splits right
#11 repts and its done

			#trim all nodes that need to be trimmed...
			for (list.index in 1:length(list.of.node.names.to.trim)) {
				#update locations/predictions of observations used during tree growth
				rows.to.update <- which(row.names(tree$frame) %in% list.of.node.names.to.trim[[list.index]])

				#create T/F vector of observations that change locations
#browser()
				observations.that.move.T.F <- names(tree$where) %in% names(list.of.trimmed.trees[[list.index]]$where)
					tree$where[observations.that.move.T.F] <- list.of.trimmed.trees[[list.index]]$where
					tree$y[observations.that.move.T.F] <- list.of.trimmed.trees[[list.index]]$y

				#update trimmed row of tree$frame
					tree$frame$n[rows.to.update] <- list.of.trimmed.trees[[list.index]]$frame$n[rows.to.update]
					tree$frame$dev[rows.to.update] <- list.of.trimmed.trees[[list.index]]$frame$dev[rows.to.update]
					tree$frame$yval[rows.to.update] <- list.of.trimmed.trees[[list.index]]$frame$yval[rows.to.update]
					tree$frame$yprob[rows.to.update,] <- list.of.trimmed.trees[[list.index]]$frame$yprob[rows.to.update,]

						#prepare for tree$frame$split
						coefs <- tree$details[[rows.to.update[1]]]$trim.sequence[dim(tree$details[[rows.to.update[1]]]$trim.sequence)[1],]
						new.oblique.split <- oblique.split.writer(	coefficients = coefs[c(TRUE,coefs[-1]!=0)],
												impurity = NULL,
												child.left.T.F = 0,
												superclass.composition = NULL)
					tree$frame$splits[rows.to.update[1],1] <- new.oblique.split$cutleft
					tree$frame$splits[rows.to.update[1],2] <- new.oblique.split$cutright

#browser()

					#tree$details and tree$details$tree.sequence
					if (is.list(new.oblique.split$details)) {
						#oblique split written as oblique split with var = ""
						tree$details[[rows.to.update[1]]]$coefficients <- new.oblique.split$details$coefficients

						#there are still oblique splits to consider
						tree$details[[rows.to.update[1]]]$trim.sequence <- tree$details[[rows.to.update[1]]]$trim.sequence[-dim(tree$details[[rows.to.update[1]]]$trim.sequence)[1],,drop=FALSE]
					} else {
						#oblique split written as axis-parallel split with var != ""
						tree$details[[rows.to.update[1]]] <- new.oblique.split$details
						tree$frame$var[rows.to.update[1]] <- new.oblique.split$variable
					}

		#		#update used entires of trim.sequence
		#		if (dim(tree$details[[rows.to.update[1]]]$trim.sequence)[1] == 1) {
		#			#set this node to "not an oblique split"
		#			tree$details[[rows.to.update[1]]]$trim.sequence <- NULL
		#		} else {
		#			#there are still oblique splits to consider
		#			tree$details[[rows.to.update[1]]]$trim.sequence <- tree$details[[rows.to.update[1]]]$trim.sequence[-dim(tree$details[[rows.to.update[1]]]$trim.sequence)[1],,drop=FALSE]
		#		}
			}

			#mark number.of.trims down by the actual number enacted
#			number.of.trims <<- number.of.trims - number.of.nodes.to.trim
			number.of.trims <<- number.of.trims - 1			#only ever trim one node at a time

		}

#look at tree after trimming or pruning
#print("just finished trimming or pruning\n");
#	subtree.leaves <- tree$frame$var == "<leaf>"
#	print(	oblique.tree:::tree.impurity(
#			yprob                                   = tree$frame$yprob[subtree.leaves,,drop=FALSE],
#			number.of.observations.at.leaves        = tree$frame$n[subtree.leaves],
#			leaf.classes                            = tree$frame$yval[subtree.leaves],
#			impurity.measure                        = "misclass")
#	)
#	print(sum(tree$y!=Pima.tr$type))
#after<-tree
#plot(after,type="u");try(text(after))
#browser()
		#housekeeping, returns nothing
		results.h[results.index] <<- min.g.of.t
		return(tree)
	}
	#######################################################################
#browser()
	#get information
	trim.impurity <- match.arg(trim.impurity)
	trim.depth <- match.arg(trim.depth)

	if (!missing(newdata)) {
		#coerce newdata into useful form
		newdata <- model.frame(	formula = as.formula(eval(tree$call$formula)),
					data = newdata)
	}

	#refit tree to training data with L1 regularization to find how variables exit the active set
	original.data <- 	model.frame(	formula = as.formula(eval(tree$call$formula)),
						data = eval(tree$call$data)
				)
	X <- tree:::tree.matrix(original.data)
	Y <- model.response(original.data)
	ylevel <- attr(tree,"ylevel")

	#everything defaulted to TRUE which denotes "an axis-parallel split has been used"
#	trim.sequence <-	lapply(	X = vector("list",length(tree$details)),
#					FUN = is.null
#				)

	#write out the changes in the active set for each oblique split
	number.of.trims <- 0				#count of the number of trims needed to reach an axis-parallel tree

	#figure out where observations fall by looking at the leaves of a rooted subtree
	leaf.node.names <- row.names(tree$frame)[tree$frame$var=="<leaf>"]

	#focus on internal nodes that apply oblique splits
#	for (details.index in which(unlist(lapply(tree$details,class))=="list")) {
#browser()
#path<-NULL
#net<-NULL
	for (details.index in which(tree$frame$var=="")) {
		#glmnet terminates before reaching saturated fits so it need not find the same fit as glm.fit() which fits to convergence even for saturated fits
		#even so, fits that stopped early still predict the training set perfectly without using all attributes in the saturated fit

		#find the leaves of subtree rooted at this node
		last.node.row.number <-	last.node.of.subtree(	tree = tree,
								subtree.root.name = as.numeric(row.names(tree$frame)[details.index])
					)
		subtree.node.names <- as.numeric(row.names(tree$frame)[details.index:last.node.row.number])
		subtree.leaf.names <- subtree.node.names[-1][subtree.node.names[-1] %in% leaf.node.names]
		subtree.leaf.row.names <- which(row.names(tree$frame) %in% subtree.leaf.names)
		observations.in.subtree.of.interest <- tree$where %in% subtree.leaf.row.names

#logistic.regression.fit
#print(details.index)
#print(c("nobs",sum(observations.in.subtree.of.interest)))
#browser()
		logistic.regression.fit <- glmnet(	x =	 	tree:::tree.matrix(	X[	observations.in.subtree.of.interest,				#observations at this node
													names(tree$details[[details.index]]$coefficients[-1]),		#attributes used at this node
													drop=FALSE
												]
									),	#model frame used at internal node during tree growth
							y =	 	as.numeric(Y[observations.in.subtree.of.interest] %in% tree$details[[details.index]]$superclass.composition),	#target variable
							family = 	"binomial",						#logistic regression
							alpha = 	1,							#no L2 penalisation wanted
							lambda.min.ratio = 	0,
#							thresh = 	1e-03,						#need 1e-04, 1e-03 is not good enough
							standardize =	FALSE)						#dont want to standardize explanatory variables to have unit variance

#		net[[details.index]] <- logistic.regression.fit
		#glmnet is perferred over glmpath as it is less particulated over lambda


#		path[[details.index]] <- glmpath(	x		=tree:::tree.matrix(	X[	observations.in.subtree.of.interest,				#observations at this node
#													names(tree$details[[details.index]]$coefficients[-1]),		#attributes used at this node
#													drop=FALSE
#												]
#										),
#							y		=as.numeric(Y[observations.in.subtree.of.interest] %in% tree$details[[details.index]]$superclass.composition),	#target variable
#							family		=binomial(),					#logistic regression
#							lambda2		=0,						#no L2 penalisation wanted
#							standardize	=FALSE,						#dont want to standardize explanatory variables to have unit variance
#							max.steps	=1000)
#}
#print(details.index)
#print(logistic.regression.fit$df)
#xxx <- FALSE									#note when df doesnt decrease monotonically
#browser()
		#we are only interested at times when active set decreases in size from what we knew before, find when these are
		active.set.of.interest <- vector("logical",length(logistic.regression.fit$df))	#all splits defaulted to ignore (FALSE)
		current.size <- logistic.regression.fit$df[length(logistic.regression.fit$df)]

		#for conformity, we dont consider the original oblique split used during tree growth. glmnet() may (or may not) stop before reaching saturated fits
		#make sure such a split is always ignored
#		if (current.size < dim(X)[2]) {
		next.size <- 1							#appropriate number of glmnet() iterations to skip
		if (current.size < length(tree$details[[details.index]]$coefficients)-1) {
			#final split isnt the same as that used during tree growth, we want it
			active.set.of.interest[length(logistic.regression.fit$df)] <- TRUE
			#next.size is effectively 1, i.e. consider the next split
		} else {
			#make sure all splits with `current.size' attributes are ignored, count the number of splits with same size
			next.size <- sum(logistic.regression.fit$df == current.size)
		}

		#look for instances when the active set decreases from what we knew before
		#the active set may increase and then decrease to a previously visited size, these "decreases" in size are not of interest, ignore them
		for (df.index in (length(logistic.regression.fit$df)-next.size):1) {
			if (logistic.regression.fit$df[df.index] < current.size) {
				#active set has reduced in size from what we knew before
				current.size <- logistic.regression.fit$df[df.index]
				active.set.of.interest[df.index] <- TRUE
			}
#if (logistic.regression.fit$df[df.index] > current.size) { xxx <- TRUE }	#note when df doesnt decrease monotonically
		}
#if (xxx) { print(details.index) }						#note when df doesnt decrease monotonically

		#an intercept only is never of interest (active set of size 1)
		active.set.of.interest[1] <- FALSE

		#collect fitted coefficients as active set changes
		tree$details[[details.index]]$trim.sequence <- 	cbind(	logistic.regression.fit$a0[active.set.of.interest],
									t(as.matrix(logistic.regression.fit$beta[,active.set.of.interest,drop=FALSE]))
								)

		#keep count of the number of trims needed to remove oblique splits
		number.of.trims <- number.of.trims + sum(active.set.of.interest)
#look at entire set of beta's found
#cbind(	logistic.regression.fit$a0,t(as.matrix(logistic.regression.fit$beta)))
	} # everything else is either a univariate split or a leaf





#browser()
#WHAT IF NOTHING TO TRIM, NEED TO WRITE THIS IN!









	#output from trim.oblique.tree() depends on what 'best' is, if it is
	#	NULL 			return the entire tree sequence
	#	a single value 		return the tree trimmed to this penalty value
	#	a vector of values 	return a summary of the trees trimmed to these values of the penalty


	if (length(best) == 1) {
		#return the trimmed tree whose penalty value is given by this value of best

		#note where the leaves are
		leaf.indicator <- tree$frame$var == "<leaf>"

		#calculate the complexity of the original tree
		current.complexity <-	oblique.tree.complexity(	tree = tree,
									subtree.internal.node.names = as.numeric(row.names(tree$frame)[!leaf.indicator])
					)

		if (current.complexity < best) {
			#cant trim the tree to the required amount
			#give a warning message
			warning("best is bigger than tree complexity")
			stop()
		} else {
			#a best tree exists, find it
			#intialise results.h to set off repeat() loop
			results.h <- 0
			results.index <- 1	#pointer for where data should be stored


			#trim just past the desired value of penalty value
			repeat {
				#calculate complexity of new tree
				leaf.indicator <- tree$frame$var == "<leaf>"
				current.complexity <-	oblique.tree.complexity(	tree = tree,
											subtree.internal.node.names = as.numeric(row.names(tree$frame)[!leaf.indicator])
							)
#print(current.complexity)
				#check if tree can actually be trimmed
				if (dim(tree$frame)[1] > 1 && current.complexity > best) {
					#the tree can be pruned and we still havent gone past the penalty value desired, find next tree in trim sequence
#browser()
					tree <- trim.once(	tree = tree,
								leaf.indicator = leaf.indicator)
				} else {
					#we either have stump or we have reached the desired tree, stop repeat() loop
					break
				}
			}
#print("here")
#browser()
			#check if we've reached the specified penalty value
			return(tree)
		}
	} else {
		#'best' is either NULL		<	>	return a summary of the entire tree sequence
		#	OR
		#a vector of 'best's 		<	>	return those trees trimmed to these penalty values

#browser()

		#create storage buckets for information from this tree sequence
		if (trim.depth == "complete") {
			#both trim and prune the tree	the maximum number of operations is 	< number.of.trims + the number of nodes
			number.of.iterations <- number.of.trims + length(tree$details)
		} else {
			#only trim the tree		the number of operations is 		= number.of.trims + 1 (the original tree)
			number.of.iterations <- number.of.trims + 1
		}
		results.complexity <- results.dev <- results.h <- vector("numeric",number.of.iterations)
#tree.save <- vector("list",number.of.iterations)

		#initialise data structures for repeat() loop
		results.h[1] <- -Inf
		results.index <- 1	#pointer for where data should be saved

		if (is.null(best)) {
			#return a summary of the entire tree sequence
			repeat {
#browser()
				#note where the leaves are
				leaf.indicator <- tree$frame$var == "<leaf>"

				#calculate the complexity of the entire tree
				results.complexity[results.index] <-	oblique.tree.complexity(	tree = tree,
													subtree.internal.node.names = as.numeric(row.names(tree$frame)[!leaf.indicator])
									)

				#calculate results.dev
				if (!missing(newdata)) {
					#newdata exists, use it to update values of $dev #and $n
					updated.tree <- predict.oblique.tree(	object =	tree,
										newdata =	newdata,
										type =		"tree",
										eps =		eps)

					if (trim.impurity == "deviance") {
						#deviance
						results.dev[results.index] <- sum(updated.tree$frame$dev[leaf.indicator])
					} else {
						#misclass
						results.dev[results.index] <- sum(updated.tree$y!=model.response(newdata))
					}
				} else {
#	browser()	tree.impurity() works fine here!
					#no newdata, find impurity of over original dataset via class.probs
					results.dev[results.index] <- tree.impurity(	yprob = tree$frame$yprob[leaf.indicator,,drop=FALSE],
											number.of.observations.at.leaves = tree$frame$n[leaf.indicator],
											leaf.classes = tree$frame$yval[leaf.indicator],
											impurity.measure = trim.impurity)
#	print(results.dev[results.index])
#	print(sum(!eval(tree$call$data)$type==tree$y))
#	the number of misclassified observations		sum(!eval(tree$call$data)$type==tree$y)
				}

#tree.save[[results.index]] <- list(tree)
#try(x11())
#try(plot(tree))
#try(text(tree))
#browser()

				#check if tree can actually be pruned
				if (trim.depth == "complete") {
					#trim and prune
					if (dim(tree$frame)[1] > 1) {
						#yes we can, find next tree in the pruning sequence
						results.index <- results.index + 1	#save tree and g.of.t to correct location
						tree <- trim.once(	tree 		= tree,
									leaf.indicator	= leaf.indicator)
					} else {
						#no, its a stump already, stop repeat() loop
						break
					}
				} else {
					#only trim
					if (number.of.trims > 0) {	#number.of.trims is altered by trim.once()
#print(paste("number.of.trims = ",number.of.trims))
						#yes we can, find next tree in the pruning sequence
						results.index <- results.index + 1	#save tree and g.of.t to correct location
						tree <- trim.once(	tree = tree,
									leaf.indicator = leaf.indicator)
					} else {
						#no, its a stump already, stop repeat() loop
						break
					}
				}
			}
#print("finished")
#browser() #inspect results
			#cut results down to size as we now know the actual number of iterations used
			results <- list(	comp = results.complexity[1:results.index],
						dev = results.dev[1:results.index],
						h = results.h[1:results.index],
						method = trim.impurity)
			attr(results,"class") <- c("trim","trim.sequence")
			return(results)
		} else {
			#initialise data structures for repeat() loop
			min.best <- min(best)

			#return a subset of the tree sequence as specified by the vector k
			repeat {
				#note where the leaves are
				leaf.indicator <- tree$frame$var == "<leaf>"

				#calculate the complexity of the entire tree
				results.complexity[results.index] <-	oblique.tree.complexity(	tree = tree,
													subtree.internal.node.names = as.numeric(row.names(tree$frame)[!leaf.indicator])
									)

				#calculate results.dev
				if (!missing(newdata)) {
					#there is newdata, use it to change the values of $dev and $n
					updated.tree <- predict.oblique.tree(	object =	tree,
										newdata =	newdata,
										type =		"tree",
										eps =		eps)

					if (trim.impurity=="deviance") {
						#deviance
						results.dev[results.index] <- sum(updated.tree$frame$dev[leaf.indicator])
					} else {
						#misclass
						results.dev[results.index] <- sum(updated.tree$y!=model.response(newdata))
					}
				} else {
					#no new data, find impurity via class.probs
					results.dev[results.index] <- tree.impurity(	yprob = tree$frame$yprob[leaf.indicator,,drop=FALSE],
											number.of.observations.at.leaves = tree$frame$n[leaf.indicator],
											leaf.classes = tree$frame$yval[leaf.indicator],
											impurity.measure = trim.impurity)
				}

				#check if we can actually prune this tree
				if (dim(tree$frame)[1] > 1 && results.complexity[results.index] > min.best) {
					#yes we can, find next tree in the pruning sequence
					results.index <- results.index + 1	#save tree and g.of.t to the correct locations
					tree <- trim.once(	tree =			tree,
								leaf.indicator =	leaf.indicator)
				} else {
					#no, its a stump already, stop this repeat() loop
					break
				}
			}

			#jig around the results and return them as asked by the vector k
			results.complexity <- results.complexity[1:results.index]
			results.dev <- results.dev[1:results.index]
			results.h <- results.h[1:results.index]

			results.complexity.star <- results.dev.star <- results.h.star <- vector("numeric",length(best))
#browser()
			for (best.index in 1:length(best)) {
				#index of interest
				correct.index <- sum(results.complexity>=best[best.index])
				results.complexity.star[best.index] <- results.complexity[correct.index]
				results.dev.star[best.index] <- results.dev[correct.index]
				results.h.star[best.index] <- results.h[correct.index]
			}

			results <- list(	comp = results.complexity.star,
						dev = results.dev.star,
						h = results.h.star,
						method = trim.impurity)
			attr(results,"class") <- c("trim","trim.sequence")
			return(results)
		}
	}
}

