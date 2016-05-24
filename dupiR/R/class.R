##Wrote by Federico Comoglio (D-BSSE, ETH Zurich) and Maurizio Rinaldi (Dipartimento di Scienza del Farmaco, University Piemonte Orientale 'A. Avogadro') 
##Last update: 06/08/2013

###########
#Class def.
###########

setClass(Class = 'Counts', representation(counts = 'integer', fractions = 'numeric', n1 = 'numeric', n2 = 'numeric', 
					  X = 'numeric', mle = 'numeric', nconst = 'numeric', posterior = 'ANY', 
					  map.p = 'numeric', map.idx = 'numeric', map = 'numeric',
					  qlow.p = 'numeric',  qlow.idx = 'integer', qlow = 'numeric', qlow.cum = 'numeric',
					  qup.p = 'numeric',  qup.idx = 'integer', qup = 'numeric', qup.cum = 'numeric',
					  gamma = 'logical'), 
		validity = function(object) {
			if( !(length(object@counts) == length(object@fractions)) ) {
				stop('the number of measurements does not match the number of fractions')
			}
			if( !all(object@fractions <= 1 & object@fractions > 0) )
				stop('Fractions not in (0,1]')
			if( any(object@counts < 0) )
				stop('Counts are negative')
			return( TRUE )
		}
)

setMethod(f = 'initialize',
		signature = 'Counts', 
		definition = function(.Object, counts, fractions) {
			if( !missing(counts) ) {
				if( sum(counts - as.integer(counts)) != 0 )
					warning('non-integer counts have been converted to integers (floor)')
			.Object@counts <-  as.integer(counts)
			if(!missing(fractions)) {
				.Object@fractions <- fractions
				.Object@X <- prod( 1 - fractions )
				mle <- round( sum(counts) / sum(fractions) )
				n1 <- round(0.5 * mle)
				n2 <- ifelse( 2 * mle == 0, round(1 / min(fractions)), 2 * mle )
				.Object@n1 <- n1
				.Object@n2 <- n2
				.Object@mle <- mle
				.Object@gamma <- FALSE
			}
			validObject(.Object)
			}
			return(.Object)
		}
)

#constructor
newCounts <- function(counts, fractions) {
	if( !missing(counts) & !missing(fractions) ) 
 		new(Class = 'Counts', counts = counts, fractions = fractions)
}
