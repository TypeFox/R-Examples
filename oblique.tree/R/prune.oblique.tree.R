`prune.oblique.tree` <- function(
	tree,						#object of class "tree"
	k = NULL,					#vector containing cost complexities to prune to, if empty, summarizes entire tree sequence, if numeric vector given, tree sequence of trees pruned to such a cost complexities is summarized and if just a single value, returns tree pruned to such a cost complexity
#	best = NULL,					#?
	newdata,					#new data with which to be used to update tree$frame
	prune.impurity = c("deviance","misclass"),	#impurity measure to use with which to compare subtrees for pruning
	penalty = c("complexity","size"),		#R_\alpha(T)=R(T)+\alpha penalty
	eps = 1e-3)
#DESCRIPTION
#	given object of class "tree", impurity measure to with which compare subtrees, will return summary of tree sequence if cost complexity k is NOT specified or it a numeric vector (greater than length 1) is specified, if of length 1, will return optimal subtree pruned to that cost complexity
#OUTPUT		impurity			numeric impurity value of subtree with such a composition of examples
{
	#######################################################################
	prune.once <- function(
		tree,						#object of class "tree" to be pruned
		leaf.indicator)		###!!! dont actually need to pass this, can just figure it out again but it wastes a little calculation
	#DESCRIPTION
	#given a tree, the next optimal tree in the optimal tree sequence is found and returned along with all g(t) for all subtrees of the tree	
	#OUTPUT		$tree=tree		object of class "tree" that resulting from pruning to the next optimal cost complexity value
	#			$g.of.t			the value of g(t) for all subtrees of the tree
	{
#x11();plot(tree);text(tree);
		#get row numbers of internal nodes 
		internal.nodes.row.numbers <- which(!leaf.indicator)

		#and the names of the leaves
		leaf.node.names <- as.numeric(row.names(tree$frame)[leaf.indicator])

		#create containers to hold information about subtrees (their leaves and g(t))
		subtree.node.names <- subtree.leaf.names <- vector("list",length(internal.nodes.row.numbers))
		g.of.t <- vector("numeric",length(internal.nodes.row.numbers))

#browser()
		#compute g.of.t
		for (subtree.index in 1:length(internal.nodes.row.numbers)) {
			subtree.row.number <- internal.nodes.row.numbers[subtree.index]
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

			#evaluate measure of misfit for subtree T.t
			R.of.T.t <- tree.impurity(	yprob = tree$frame$yprob[subtree.leaves.rows.T.F,,drop=FALSE],
							number.of.observations.at.leaves = tree$frame$n[subtree.leaves.rows.T.F],
							leaf.classes = tree$frame$yval[subtree.leaves.rows.T.F],
							impurity.measure = prune.impurity)

			#evaluate measure of misfit for subtree T.t
			if (penalty == "size") {
				#...when cost-complexity penalises number of leaves
				penalty.of.R.of.T.t <- sum(subtree.leaves.rows.T.F)
			} else {
				#...when cost-complexity penalises complexity of branches leading to each leaf
				penalty.of.R.of.T.t <- 	oblique.tree.complexity(	tree = tree,
											subtree.internal.node.names = subtree.node.names[[subtree.index]][!(subtree.node.names[[subtree.index]] %in% leaf.node.names)]
							)
			}

			#evaluate measure of misfit for pruned node t
			R.of.t <- tree.impurity(	yprob = tree$frame$yprob[subtree.row.number,,drop=FALSE],
							number.of.observations.at.leaves = tree$frame$n[subtree.row.number],
							leaf.classes = tree$frame$yval[subtree.row.number],
							impurity.measure = prune.impurity)

			#evaluate g.of.t (penalty.of.R.of.t is 1 as we prune to a node)
			g.of.t[subtree.index] <- (R.of.t - R.of.T.t) / (penalty.of.R.of.T.t - 1)
		}
#print(g.of.t)
#browser()

		#nodes with min(g.of.t) need to be pruned
		min.g.of.t <- min(g.of.t)
#print(which.min(g.of.t))

		#prune nodes that need to be pruned
		nodes.to.prune.T.F <- g.of.t == min.g.of.t
		if (sum(nodes.to.prune.T.F) > 0) {
			#node(s) needs pruning, alter	$frame and row.names($frame) to reflect the updated tree, 
			#				$where to the new locations using the old node numbering system
			#				$y with the correct new predictions
			#whilst keep a running idea of which nodes have been pruned to allow removal of unused entries in mat and $details in one go
			tree <- oblique.tree.prune.nodes(	tree = tree,
								list.of.node.names.to.prune = subtree.node.names[nodes.to.prune.T.F]
				)
		}

		#housekeeping, returns nothing
		results.k[results.index] <<- min.g.of.t
		return(tree)
	}
	#######################################################################

	#get information
	prune.impurity <- match.arg(prune.impurity)
	penalty <- match.arg(penalty)

	if (!missing(newdata)) {
		#coerce newdata into useful form
		newdata <- model.frame(	formula = as.formula(eval(tree$call$formula)),
					data = newdata)
	}

	#output from prune.oblique.tree() k may be NULL (just give me the tree sequence), a single value (give me a pruned tree with this value of k) or a vector of values (give me a summary of the trees pruned to these values of k), do the appropriate thing
	if (length(k)==1) {
		#return the rooted subtree pruned to cost complexity k (which here is just a number)
		#intialise results.k to set off repeat() loop
		results.k <- vector("numeric",dim(tree$frame)[1])
		results.k[1] <- -Inf
		results.index <- 1

		#prune to the desired cost complexity
		old.tree <- tree					#in the case where pruning to -Inf
		repeat {
			#check if we can actually prune this tree
			if (dim(tree$frame)[1] > 1 && k != -Inf && round(results.k[results.index]-k,6) <= 0) {
				#yes we can, find next tree in the pruning sequence
				results.index <- results.index + 1	#save tree and g.of.t to the correct locations
				old.tree <- tree			#save the old tree and look ahead with 'tree'
				tree <- prune.once(	tree = tree,
							leaf.indicator = tree$frame$var=="<leaf>")
			} else {
				#no, we've either hit a stump or its a stump already, stop this repeat() loop
				break
			}
		}

		#check if we've reached the specified k yet even though we have a root tree? put best tree in 'tree'
		if (round(results.k[results.index]-k,6) <= 0 || k==-Inf) {
			return(tree)
		} else {
			return(old.tree)
		}
	} else {
		#k is either NULL	<	>	return a summary of the tree sequence
		#	OR
		#a vector of k's 	<	>	return those trees denoted by this cost-complexity sequence 

		#create data stores for information about tree sequence
		leaf.indicator <- tree$frame$var=="<leaf>"			#make a note of where the leaves are in tree$frame
		container.size <- sum(leaf.indicator)				#upper limit on number of saved rooted-subtrees required to get to tree stump
		results.complexity <- results.dev <- results.k <- vector("numeric",container.size)

		#initialise data structures for repeat() loop
		results.k[1] <- -Inf					#save results.k
		results.index <- 1

		if (is.null(k)) {
			#return a summary of the entire tree sequence
			repeat {
				#save results.complexity
				leaf.indicator <- tree$frame$var=="<leaf>"			#note where leaves are in tree$frame
				results.complexity[results.index] <- sum(leaf.indicator)

				#save results.dev
				if (!missing(newdata)) {
					#newdata exists, use it to update values of $dev #and $n
					updated.tree <- predict.oblique.tree(	object = tree,
										newdata = newdata,
										type = "tree",
										eps = eps)
					if (prune.impurity=="deviance") {
						#deviance
						results.dev[results.index] <- sum(updated.tree$frame$dev[leaf.indicator])
					} else {
						#misclass
						results.dev[results.index] <- sum(updated.tree$y!=model.response(newdata))
					}
				} else {
					#no newdata, find impurity of over original dataset via class.probs
					results.dev[results.index] <- tree.impurity(	yprob = tree$frame$yprob[leaf.indicator,,drop=FALSE],
											number.of.observations.at.leaves = tree$frame$n[leaf.indicator],
											leaf.classes = tree$frame$yval[leaf.indicator],
											impurity.measure = prune.impurity)
				}

				#check if tree can actually be pruned
				if (dim(tree$frame)[1]>1) {
					#yes we can, find next tree in the pruning sequence
					results.index <- results.index + 1	#save tree and g.of.t to the correct locations
					tree <- prune.once(	tree=tree,
								leaf.indicator=leaf.indicator)
				} else {
					#no, its a stump already, stop repeat() loop
					break
				}
			}

			#cut results.* down to size as we know the actual number of pruning operations required to reach tree stump
			results <- list(	size=results.complexity[1:results.index],
						dev=results.dev[1:results.index],
						k=results.k[1:results.index],
						method=prune.impurity)
			attr(results,"class") <- c("prune","tree.sequence")
			return(results)
		} else {
			#initialise data structures for repeat() loop
			max.k <- max(k)

			#return a subset of the tree sequence as specified by the vector k
			repeat {
				#save results.complexity
				leaf.indicator <- tree$frame$var=="<leaf>"			#make a note of where the leaves are in tree$frame
				results.complexity[results.index] <- sum(leaf.indicator)

				#save results.dev
				if (!missing(newdata)) {
					#there is newdata, use it to change the values of $dev and $n
					updated.tree <- predict.oblique.tree(	object=tree,
										newdata=newdata,
										type="tree",
										eps=eps)
					if (prune.impurity=="deviance") {
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
											impurity.measure = prune.impurity)
				}

				#check if we can actually prune this tree
				if (dim(tree$frame)[1] > 1 && results.k[results.index] < max.k) {
					#yes we can, find next tree in the pruning sequence
					results.index <- results.index + 1	#save tree and g.of.t to the correct locations
					tree <- prune.once(	tree=tree,
								leaf.indicator=leaf.indicator)
				} else {
					#no, its a stump already, stop this repeat() loop
					break
				}
			}

			#jig around the results and return them as asked by the vector k
			results.complexity <- results.complexity[1:results.index]
			results.dev <- results.dev[1:results.index]
			results.k <- results.k[1:results.index]

			results.complexity.star <- results.dev.star <- results.k.star <- vector("numeric",length(k))

			for (k.index in 1:length(k)) {
				#index of interest
				correct.k <- sum(results.k <= k[k.index])
				results.complexity.star[k.index] <- results.complexity[sum(results.k <= k[k.index])]
				results.dev.star[k.index] <- results.dev[sum(results.k <= k[k.index])]
				results.k.star[k.index] <- k[k.index]
			}

			results <- list(	size=results.complexity.star,
						dev=results.dev.star,
						k=results.k.star,
						method=prune.impurity)
			attr(results,"class") <- c("prune","tree.sequence")
			return(results)
		}
	}
}

