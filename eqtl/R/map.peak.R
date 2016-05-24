#####################################################################
#
# map.peak.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/eqtl package
# Contains: map.peak 
#
######################################################################

######################################################################
#
# map.peak: Resume maximum LOD peak position from peak object.
#
######################################################################

`map.peak` <-
function(peak)
{

if ( !all(attr(peak,'class',exact=TRUE) %in% c('peak','list')) )
	stop("Input should have class \"peak\".")

trait <- NA
chr <- NA
peak.cM <- NA

for (i in seq(names(peak))){
	for ( j in 1:length(peak[[i]])){
		if ( is.data.frame(peak[[i]][[j]]) ) {
			for ( z in seq(length(as.vector(peak[[i]][[j]]$mname.peak))) ){
				trait <- c(trait,names(peak[i]))
				chr <- c(chr,names(peak[[i]][j]))
				peak.cM <- c(peak.cM,peak[[i]][[j]]$peak.cM[z])
			}
		}
	}
}

res <- data.frame(trait[-1],chr[-1],peak.cM[-1])
attributes(res)$names <- c('trait','chr','cM')
return(res);

}

