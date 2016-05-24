
############################
# helper function significance level
label_significance_level <- function( values , levels , labels ){
		ix <- sort( levels, index.return=TRUE)$ix
		levels <- levels[ix]
		labels <- labels[ix]
		NL <- length(levels)
		l1 <- ""
		values[ is.na(values) ] <- 1.2
		for (ll in 1:NL){
			l1 <- ifelse( values < levels[NL-ll+1] , labels[NL-ll+1] , l1)
						}
		return(l1)
				}
###########################################				