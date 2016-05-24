#####################################################################
#
# wash.covar.R
#
# copyright (c) 2008-3, Ahmid A Khalili
# 
# last modified Jul, 2008
# first written Mar, 2008
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/eqtl package
# Contains: wash.covar
#
######################################################################

######################################################################
#
# wash.covar: Erase additive covariates LOD peaks on LOD curve.
#
######################################################################

`wash.covar` <-
function( scanone, covar, window.size=20)
{

	if ( !all(attr(scanone,'class',exact=TRUE) %in% c('scanone','data.frame')) )
        	stop("Input should have class \"scanone\".")

	if ( !is.data.frame(covar) & any(!(names(covar) %in% c('trait','chr','cM')))  )
        	stop("Argument covar misspecified: covar should be a data.frame with columns names 'trait','chr','cM'")
	if ( any(!(levels(covar$chr) %in% levels(scanone$chr))) )
        	stop("Argument covar misspecified: wrong covar$chr value")
	if ( !is.vector(covar$cM) & !is.numeric(covar$cM) & any(covar$cM < 0) )
        	stop("Argument covar misspecified: wrong covar$cM value")
	if ( length(levels(covar$trait)) > ncol(scanone)-2 ) # as.vector ???
        	stop("Argument covar misspecified: wrong covar$trait value")
	if ( !is.vector(window.size) & !is.numeric(window.size) & length(window.size)!=1 & window.size<0 )
        	stop("Argument misspecified: wrong window.size value (centiMorgan)")

	chr <- as.numeric(scanone$chr)
	pos <- as.numeric(scanone$pos)
	trait <- names(scanone)

	final_scanone <- scanone[1:2]
	class(final_scanone) <- 'data.frame'

	if ( length(trait)>3 ){
		for ( i in seq(length(covar$trait))){

			# INFORMATIVE SCREEN MESSAGE
			# cat("Searching trait",as.vector(covar$trait[i]),": ")

			z <- 3

			while ( z <= length(trait)){
				if ( toupper(covar$trait[i]) == toupper(trait[z]) ){

					# INFORMATIVE SCREEN MESSAGE 
					# cat(trait[z],"found z:",z,"\n")

					bool <- (pos >= covar$cM[i]-window.size) & (pos <= covar$cM[i]+window.size) & (chr == covar$chr[i])
					if ( !any(bool) ) stop("ERROR in boolean line 34")

					if ( i == 1 ){
						new_scanone <- scanone[z]
						new_scanone[bool,] <- rep(0,length(bool[bool==TRUE]))
						class(new_scanone) <- 'data.frame'
						final_scanone <- cbind(final_scanone,new_scanone)
					} else {
						if ( as.vector(covar$trait[i]) == as.vector(covar$trait[i-1]) ){
							final_scanone[bool,length(names(final_scanone))] <- rep(0,length(bool[bool==TRUE]))
						} else {
							new_scanone <- scanone[z]
							new_scanone[bool,] <- rep(0,length(bool[bool==TRUE]))
							class(new_scanone) <- 'data.frame'
							final_scanone <- cbind(final_scanone,new_scanone)
						}
					}

					if ( trait[z] != trait[z+1] & z<length(trait)) break	
				}
				z <- z+1
				next
			}
		}
	} else {

		# INFORMATIVE SCREEN MESSAGE
		# cat("WARNING: only one trait")

		for ( i in seq(length(covar$trait))){
		
			bool <- (pos >= covar$cM[i]-window.size) & (pos <= covar$cM[i]+window.size) & (chr == covar$chr[i])

			if ( i == 1 ){
				new_scanone <- scanone[3]
				new_scanone[bool,] <- rep(0,length(bool[bool==TRUE]))
				class(new_scanone) <- 'data.frame'
				final_scanone <- cbind(final_scanone,new_scanone)
			} else {
				if ( as.vector(covar$trait[i]) == as.vector(covar$trait[i-1]) ){
					final_scanone[bool,3] <- rep(0,length(bool[bool==TRUE]))
				} else {
					new_scanone <- scanone[3]
					new_scanone[bool,] <- rep(0,length(bool[bool==TRUE]))
					class(new_scanone) <- 'data.frame'
					final_scanone <- cbind(final_scanone,new_scanone)
				}
			}
		}
	}

	attributes(final_scanone) <- attributes(scanone)
	return(final_scanone)
}

