#####################################################################
#
# drop.peakfeat.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/eqtl package
# Contains: drop.peakfeat
#
######################################################################

######################################################################
#
# drop.peakfeat: Erase chosen peak features informations from
#                       a \code{peak} object.
#
######################################################################

`drop.peakfeat` <-
function(peak,feat)
{

	require(qtl)

	if ( !all(attr(peak,'class',exact=TRUE) %in% c('peak','list')) )
		stop("Input should have class \"peak\".")
	if ( missing(peak) )
		stop("Argument 'peak' unspecified.")
	if ( !is.vector(feat) || !any(feat %in% attr(peak,'features',exact=TRUE)) )
		stop("Argument 'peak' misspecified: Expecting a vector containing the 'peak' features attributes (",attr(peak,'features',exact=TRUE),")")
	for (i in 1:length(peak)){ 
		for (y in 1:length(peak[[i]])){
			if (!is.na(peak[[i]][y])){
				for ( f in feat ){
					if ( f %in% names(peak[[i]][[y]] ) ){
						col <- grep(f,names(peak[[i]][[y]]))
						peak[[i]][[y]] <- peak[[i]][[y]][-col]
					} else print("feat is not defined in peak")
				}
			}
		}
	}

	attributes(peak)$features <-  attributes(peak)$features[ ! attr(peak,'features') %in% feat ]
	try(return(peak),silent=FALSE)
}

