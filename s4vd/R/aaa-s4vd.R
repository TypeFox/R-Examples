setClass('BCs4vd',
		contains = 'BiclustMethod',
		prototype = prototype(
				biclustFunction = function(X,...){s4vd(X,...)}))



BCs4vd <- function() {
	return(new('BCs4vd'))
}


setClass('BCssvd',
		contains = 'BiclustMethod',
		prototype = prototype(
				biclustFunction = function(X,...){ssvdBC(X,...)}))

BCssvd <- function() {
	return(new('BCssvd'))
}

