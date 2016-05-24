##Script generated in:
# 2011
# 10:51:13 PM
#by: 
# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################


reduceStructure <- function(knot, algorithm = 'AB') {
			stopifnot( is(knot, 'Knot') )
			if( missing(algorithm) || !( algorithm %in% c('AB', 'MSR')) )
				return('algorithm is missing or non correct: type \'AB\' for AlexanderBriggs or \'MSR\' for Minimal Structure Reduction')
			if(algorithm == 'AB'){
				reduced <- AlexanderBriggs(knot@points3D, knot@ends)
				knot@points3D <- reduced$points3D
				knot@ends <- reduced$ends
				return(knot)
			}
			if(algorithm == 'MSR') {
				reduced <- msr(knot@points3D, knot@ends)
				knot@points3D <- reduced$points3D
				knot@ends <- reduced$ends
				return(knot)
			}
		}

computeInvariant <- function(knot, invariant, ...) {
	stopifnot( is(knot, 'Knot') )
	implemented <- c('Alexander', 'Jones', 'HOMFLY', 'LK')
	if( missing(invariant) || !( invariant %in% implemented) )
				return('invariant not specified')
	switch(invariant,
		'Alexander' = {
						reduced <- reduceStructure(knot, algorithm = 'AB')
						inv <- mVA(reduced@points3D, reduced@ends, ...)
						return(inv) 
					   },
		'Jones' = {
						reduced <- reduceStructure(knot, algorithm = 'AB')
						tmp <- skeinIterator(reduced@points3D, reduced@ends)
						inv <- HOMFLYpolynomial(tmp$leaves, tmp$tree, ...)
						inv <- HOMFLY2Jones( inv )
						return(inv)
					},
		'HOMFLY' = {
						reduced <- reduceStructure(knot, algorithm = 'AB')
						tmp <- skeinIterator(reduced@points3D, reduced@ends)
						inv <- HOMFLYpolynomial(tmp$leaves, tmp$tree, ...)
						return(inv)
					},
		'LK' = {
						if(!identical( knot@ends, numeric(0) )) {
							reduced <- reduceStructure(knot, algorithm = 'AB')
							inv <- linkingNumber(knot@points3D, knot@ends)
							return(inv)
						}
						else return('Linking number computation requires a link as a input. ends are missing.')
					}
			)
}

closeAndProject <- function(protein, ...) {
	stopifnot( is(protein, 'Knot') )
	tmp <- centroidClosure(protein@points3D, ...) #can specify w
	protein@points3D <- PCAProjection(tmp)
	return(protein)
}

plot3D <- function(knot, ...) {
	plotKnot3D(knot@points3D, knot@ends, ...)
}


