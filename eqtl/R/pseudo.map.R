#####################################################################
#
# pseudo.map.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/eqtl package
# Contains: pseudo.map
#
######################################################################

######################################################################
#
# pseudo.map: List the markers and pseudo-markers genetic positions
#
######################################################################

`pseudo.map` <-
function(cross)
{

	require(qtl)

	if ( class(cross)[2] != "cross" ) stop("peak should have class \"cross\"\n")
	if ( ! 'prob' %in% names(cross$geno[[1]]) ) stop('First running calc.genoprob')

	map <- NA
	for ( i in seq(names(cross$geno)) ) map <- c(map, attr(cross$geno[[i]]$prob,'map') )

	return(map[-1])

}

