setClass("octav", representation("data.frame"), validity = function(object) {
		 if(dim(object)[2] != 3) return("octav objects must have 3 columns")
		 if(! (is.numeric(object[,1]) & is.numeric(object[,2]) & is.numeric(object[,3]))) return ("All columns must be of numeric or integer class")
		 if(any(object[,1]!=object[order(object[,1]),1]) | any(object[,2]!=object[order(object[,2]),2])) return("The octav and upper columns of the octav class must be ordered")
		 if(any(object[,2] < 0) | any(object[,3] < 0)) return("No negative values are allowed in the upper and Freq columns")
		 if(any(abs(object[,2] - 2^object[,1]) > 1e-8)) return("The upper column must be made of powers of the octav column")
		 TRUE
}
)

setClass("rad", representation("data.frame"), validity = function(object) {
		 if(dim(object)[2] != 2) return("rad objects must have 2 columns")
		 if(! (is.numeric(object[,1]) & is.numeric(object[,2]))) return ("All columns must be of numeric or integer class")
		 ab = object[,2]
		 if(any(is.na(ab))) return("NAs generated in the rad object!!")
		 if(any(ab!=ab[order(-ab)])) return("The abund column of the rad class must be in descending order")
		 if(any(object[,1] < 0) | any(ab < 0)) return("No negative values are allowed in the rad object")
		 TRUE
}
)

setClass("fitsad", representation("mle2", sad="character", distr="character", trunc="numeric"))
setClass("fitrad", representation("mle2", rad="character", distr="character", trunc="numeric", rad.tab="rad"))

distr.depr <- "The 'distr' slot of fitrad and fitsad objects have been deprecated. Please see ?distr"

#' Summary for fitsad/fitrad calls
#' 
#' This function works almost exactly as bbmle's summary.mle2, but it includes a "fixed parameters" 
#' line for models with fixed parameters, such as \code{\link{fitls}} or \code{fitvolkov}.
#' NOTICE that the summary.mle2 is redefined in this package, as this class is not exported by bbmle.
#' @rdname summary.sads-class
setClass("summary.mle2", representation(call = "language", coef = "matrix",m2logL = "numeric", fixed="numeric"))
#' @rdname summary.sads-class
setClass("summary.sads", representation("summary.mle2", fixed="numeric"))
