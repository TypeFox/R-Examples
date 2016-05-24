`oblique.split.writer` <- function(			#	given
	coefficients,					#		a) a vector of named coefficients of the oblique split of interest
	impurity,					#		b) the impurity of the split specified by the vector of coefficients
	child.left.T.F,					#		c) a vector of TRUEs and FALSEs of observations partitioned by this split into child nodes
	superclass.composition)				#		d) the composition of the superclass that led to this split so that the targets of fitted model can be generated
{							#	oblique splits are written nicely when p>=2 so they can be easily handled
	#return a list
	split <- NULL

	#oblique split can be considered as axis parallel splits and so need to cater for it here
	split$impurity <- impurity
	if (length(coefficients) == 2) {
		#have univariate split
		cutoff.point <- round(-coefficients[1]/coefficients[2],6)
		split$cutleft <- paste("<",cutoff.point,sep="")
		split$cutright <- paste(">",cutoff.point,sep="")
		split$details <- as.numeric(cutoff.point)

		#axis-parallel splits masquerading as oblique splits need to be handled with care, child.left.T.F is c+aX<0 is not X<-c/a, it depends on sign of 'a'
		#the linear predictor of oblique splits has evaluated c+aX<0 to give child.left.T.F
		if (coefficients[2] < 0) {
			#when a<0, we have X>-c/a, which is _not_ what we can, need to invert the allocations to child nodes
			split$child.left.T.F <- !child.left.T.F
		} else {
			#when a>0, we have X<-c/a, which is what we want
			split$child.left.T.F <- child.left.T.F
		}
		split$variable <- names(coefficients)[2]
	} else {
		#have linear-multivariate split
		#make cutoff.point look nice, no -+ and spaces
		cutoff.point.sign <- rep("+",(length(coefficients)-1))	#default to +
		cutoff.point.sign[coefficients[-1]<0] <- "-"		#change correct terms over to - to look nice
		cutoff.point <- vector("character",3*length(coefficients)-2)
		cutoff.point[1] <- round(coefficients[1],2)		#will change using "for" loop
		for (cutoff.point.index in 1:(length(coefficients)-1)) {
			cutoff.point[3 * cutoff.point.index-1] <- cutoff.point.sign[cutoff.point.index]
			cutoff.point[3 * cutoff.point.index] <- round(abs(coefficients[1+cutoff.point.index]),2)
			cutoff.point[3 * cutoff.point.index+1] <- names(coefficients[1+cutoff.point.index])
		}
		cutoff.point <- paste(cutoff.point,collapse="")
		split$cutleft <- paste(cutoff.point,"<0",sep="")
		split$cutright <- paste(cutoff.point,">=0",sep="")
		split$details$coefficients <- coefficients
		split$details$superclass.composition <- superclass.composition			#store composition of left superclass for use in prune.oblique.tree()
		split$child.left.T.F <- child.left.T.F
		split$variable <- ""
	}
	return(split)
	#OUTPUT		split	$impurity
	#			$child.left.T.F
	#			$details	$coefficients
	#					$superclass.composition
	#			$variable
	#			$cutleft
	#			$cutright
}

