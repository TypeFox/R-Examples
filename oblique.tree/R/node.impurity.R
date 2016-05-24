`node.impurity` <- function(			#given
	class.probabilities,			#	a) a table of class probabilities at a node
	impurity.measure = c(	"deviance",	#	b) impurity measure to use when growing tree
				"gini"))	#the impurity of a node is returned
{
	#evaluate impurity when given class probabilities
	if (impurity.measure == "deviance") {
		#deviance
		class.probabilities[class.probabilities == 0] <- 1			#change any class.probs from 0 to 1 to make 0*log(0)=0 by default
		impurity <- -2*	sum(	class.probabilities * 
					log(class.probabilities)
				)
	} else if (impurity.measure == "gini") {
		#gini
		impurity <- 1 - sum(class.probabilities^2)
	}

	return(impurity)
}

