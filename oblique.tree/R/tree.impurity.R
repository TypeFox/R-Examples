`tree.impurity` <- function(
	yprob,						#table containing probability of classes of all leaves in subtree
	number.of.observations.at.leaves,		#vector containing the number of examples of all leaves in subtree
	leaf.classes,
	impurity.measure = c("deviance","misclass"))	#impurity measure to use with which to compare subtrees
#DESCRIPTION
#	given class probability values, frequency, predicted class and node impurity measure, all leaves in subtree are considered to return an impurity measure with which to compare subtrees	
#OUTPUT		impurity		numeric impurity value of subtree with such a composition of examples
{
	if (impurity.measure == "deviance") {
		#R(T) measured by deviance of observations
		yprob[yprob==0] <- 1
		impurity <- -2 * sum(	yprob * log(yprob) * number.of.observations.at.leaves)
	} else {
#browser()
		#R(T) measured by misclassified observations
		impurity <- 0
		for (leaf.index in 1:length(number.of.observations.at.leaves)) {
			#sum over all leaves to find number of misclassified examples with respect to y
			impurity <- impurity + 	sum(	yprob[leaf.index,-as.numeric(leaf.classes[leaf.index])]
						) * number.of.observations.at.leaves[leaf.index]
#print(c(leaf.index,impurity))
		}
	}
	return(impurity)
}

