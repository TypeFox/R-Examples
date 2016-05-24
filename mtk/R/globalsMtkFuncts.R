# Mexico Toolkit
#
# version   : 0.03
# date		: Sept 2011
# MAJ   	: 
# licence	: GPL

# Author(s) : Hervé Richard INRA MIA BioSP
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev::                   $ : revision number of the last spread
#' $Author::                $ : author of the last spread
#' $Date::                  $ : date of the last spread


#=======================================================================
# MISCELLANEOUS FUNCTIONS
#=======================================================================

#'  Calculates the quantiles of a univariate distribution
#'  @title The quantiles of a univariate distribution
#'  @param pvalues a vector of probability values
#'  @param distribName a character string: name of a probability distribution
#'  @param distribParameters a list of parameters of the distribution
#'  @param shrink a scalar eqn{<=1} to determine how to shrink the pvalues - used when
#'         the quantiles are infinite for pvalues equal to 0 or 1
#'  @examples
#'   Quantiles(seq(0,1,length=11),"unif",list(min=8,max=10))
#'   Quantiles(seq(0,1,length=11),"unif",list(min=8,max=10),shrink=0.5)
#'   Quantiles(seq(0,1,length=11),"norm",list(mean=0, sd=1),shrink=0.5)
#'  @export Quantiles

Quantiles <- function(pvalues, distribName, distribParameters, shrink=0.95)
{
	quantileFunctionName <- paste("q",distribName,sep="")
	qvalues <- eval( do.call(quantileFunctionName,
					c(list(p=pvalues),distribParameters))
	)
	
	# warning if there are infinite quantiles
	if( (Inf %in% abs(qvalues)) ){
		# shrinkage
		message <- paste("Infinite values in the design.",
				"P-values shrinked by",
				shrink
		)
		warning(message)
		pvalues <- 0.5 + (pvalues-0.5)*shrink
		qvalues <- eval( do.call(quantileFunctionName,
						c(list(p=pvalues),distribParameters))
		)
	}
	# output
	return(qvalues)
}



