# 
# Ingmar Viser, 23-3-2008
# 

setClass("llratio",
	representation(
		value="numeric",
		df="numeric"
	)
)

llratio <- function(basemodel,constrainedmodel,...) {
	llbase <- logLik(basemodel)
	llcon <- logLik(constrainedmodel)
	llr <- 2*(llbase-llcon)
	df <- attributes(llbase)$df-attributes(llcon)$df
	return(new("llratio",value=as.numeric(llr),df=df))
}

setMethod("show","llratio",
	function(object) {
		pvalue=round(pchisq(object@value,object@df,lower.tail=FALSE),3)
		cat("log Likelihood ratio (chi^2): ", round(object@value,3), " (df=",object@df,"), p=",pvalue,".\n",sep="")
	}
)
