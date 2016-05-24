#' Conditional Element Selection for \code{ff} vectors.
#'
#' Similar as \code{ifelse} in the base package but only works with yes and no as \code{ff} vectors.
#'
#' @export
#' @example ../examples/ffifelse.R
#' @param test logical or boolean \code{ff} vector
#' @param yes an \code{ff} vector with return values for true elements of test. If too short, their elements are recycled.
#' @param no an \code{ff} vector with return values for false elements of test. If too short, their elements are recycled.
#' @return An ff vector of the same length as \code{test}. 
#' @seealso \code{\link[base]{ifelse}}
ffifelse <- function(test, yes, no){
	if(!vmode(test) %in% c("boolean","logical") || !is.ff(test)){
		stop("test needs to be a logical/boolean ff vector")
	}
	if(length(yes) != length(test) || !is.ff(yes)){
		yes <- ff(yes, length=length(test))
	}
	if(length(no) != length(test) || !is.ff(no)){
		no <- ff(no, length=length(test))
	}
	dat <- ffdf(fftest = test, ffyes = yes, ffno = no)
	if(is.factor.ff(dat$ffyes) && is.factor.ff(dat$ffno)){
		result <- with(dat, ifelse(fftest, as.character(ffyes), as.character(ffno)))
	}else if(is.factor.ff(dat$ffyes) && !is.factor.ff(dat$ffno)){
		result <- with(dat, ifelse(fftest, as.character(ffyes), ffno))	
	}else if(!is.factor.ff(dat$ffyes) && is.factor.ff(dat$ffno)){
		result <- with(dat, ifelse(fftest, ffyes, as.character(ffno)))
	}else{
		result <- with(dat, ifelse(fftest, ffyes, ffno))	
	}
	result	
}

