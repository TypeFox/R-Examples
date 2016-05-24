#####################################################################
#
# peak.2.array.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/eqtl package
# Contains: peak.2.array 
#
######################################################################

######################################################################
#
# peak.2.array: Build QTL results array from peak object
#
######################################################################

`peak.2.array` <-
function(peak)
{

	if ( !all(attr(peak,'class',exact=TRUE) %in% c('peak','list')) )
		stop("Input should have class \"peak\".")

	for ( i in seq(names(peak)) ){
		for (j in 1:length(peak[[i]])) {
			if ( is.data.frame(peak[[i]][[j]])){
				feat <- names(peak[[i]][[j]])
				break
			}
		}
		if (is.data.frame(peak[[i]][[j]])) break
	}

	cat('peaks features:\n')
	print(feat)
	cat('\n')

	array <- data.frame(t(rep(NA,length(feat)+2)))

	for ( i in seq(names(peak)) ){
		for ( j in 1:length(peak[[i]]) ) {

			if ( is.data.frame(peak[[i]][[j]]) ) {
				for ( z in seq(length(as.vector(peak[[i]][[j]]$mname.peak))) ){
					w <- NA
					for ( y in 1:length(feat) ) w <- c(w,as.vector(peak[[i]][[j]][[y]][z]))
					array <- rbind(	array,data.frame(t(c(	names(peak[i]),
										names(peak[[i]][j]),
										w[-1]
									)))
							)
				}
			} else {
				array <- rbind(	array,data.frame(t(c(	names(peak[i]),
									names(peak[[i]][j]),
									rep(NA,length(feat))
								)))
						)
			}
		}
	}

	attributes(array)$names <- c('trait','chr',feat)
	array <- array[-1,]
	attributes(array)$class <- c('peak.array','data.frame')
	return(array)
}

