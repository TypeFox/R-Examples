##Script generated in:
# 2011
# 9:55:08 PM
#by: 
# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

#makeExampleKnot (k governs the output. k = TRUE returns a knot, FALSE returns a link.
makeExampleKnot <- function(k = TRUE) {
	if(k) {
		#Rolfsen.table <- NULL; rm(Rolfsen.table)
		#data(Rolfsen.table, package='Rknots')
		knot <- Rknots::Rolfsen.table[[ sample(1 : length(Rknots::Rolfsen.table), size = 1) ]]
		return( knot )
	}
	else {
		#link.table <- NULL; rm(link.table)
		#data(link.table, package='Rknots')
		link <- Rknots::link.table[[ sample(1 : length(Rknots::link.table), size = 1) ]]
		return( link )
	}
}

makeExampleProtein <- function() {
	fn <- system.file( "extdata", "2K0A.pdb", package = "Rknots")
	protein <- loadProtein( fn )
	return(protein)
}
