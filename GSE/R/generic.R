## set a subclass for all the classes created
setClass("CovRobMiss", representation(mu="vector", 
			S="matrix", 
			call="language", 
			estimator="character", 
			x="matrix",
			pmd="vector", 
			pmd.adj="vector",
			p="numeric",
			pu = "vector")) 
setClass("CovRobMissSc", representation(
			sc="numeric"), contains="CovRobMiss")				
setClass("emve", representation(), contains="CovRobMissSc")
setClass("HuberPairwise", representation(R="matrix"), contains="CovRobMiss")
setClass("GSE", representation(
			mu0 = "vector",
			S0 = "matrix",
			weights = "vector",
			weightsp = "vector",
			ximp = "matrix",
			iter = "numeric",
			eps="numeric"), contains="CovRobMissSc")
setClass("SummaryCov", representation(obj="CovRobMiss", evals="list"))
setClass("TSGS", representation(xf="matrix"), contains="GSE")
setGeneric("getFiltDat", function(object) standardGeneric("getFiltDat"))
setMethod("getFiltDat", "TSGS", function(object) object@xf)


## S4 method of printing 
setMethod("show", "CovRobMiss", function(object){
	digits = max(3, getOption("digits") - 3)
    cat("\nCall:\n", deparse(object@call), "\n\n", sep = "")
    cat("-> Estimator: ", object@estimator, "\n")
    cat("\nEstimate of Location: \n")
    print.default(format(object@mu, digits = digits), print.gap = 2, quote = FALSE)
    cat("\nEstimate of Covariance: \n")
    print.default(format(object@S, digits = digits), print.gap = 2, quote = FALSE)
    invisible(object)
})

## S4 method of summary
## in addition to the basics as in printing, also include:
## Proportion of missingness, Any rows with completely missing data,
## Eigenvalues of estimated covariance,
## Partial mahalanobis distances
setMethod("summary", "CovRobMiss", function(object){
	## obtain eigenvalues 
	new("SummaryCov", obj=object, evals=eigen(object@S))
})

setMethod("show", "SummaryCov", function(object){
	digits = max(3, getOption("digits") - 3)
    cat("\nCall:\n", deparse(object@obj@call), "\n\n", sep = "")
    cat("-> Estimator: ", object@obj@estimator, "\n")
    cat("\nEstimate of Location: \n")
    print.default(format(object@obj@mu, digits = digits), print.gap = 2, quote = FALSE)
    cat("\nEstimate of Covariance: \n")
    print.default(format(object@obj@S, digits = digits), print.gap = 2, quote = FALSE)
	cat("\nProportion of missingness: \n")
	print.default(format( mean(is.na(object@obj@x)), digits = digits), print.gap = 2, quote=FALSE)
	cat("\nRow(s) with completely missing data: \n")
	print.default(format( as.numeric(which( rowSums(is.na(object@obj@x)) == ncol(object@obj@x))), digits = digits), print.gap = 2, quote=FALSE)
	cat("\nEigenvalues of covariance matrix: \n")
    print.default(format(object@evals$values, digits = digits), print.gap = 2, quote = FALSE)
    cat("\nAdjusted squared partial Mahalanobis Distances: \n")
    print.default(format(as.vector(object@obj@pmd.adj), digits = digits), print.gap = 2, quote = FALSE, max=10)
    invisible(object)
})

setGeneric("getDistAdj", function(object) standardGeneric("getDistAdj"))
setMethod("getDistAdj", "CovRobMiss", function(object){
	pmd.adj <- object@pmd.adj
	## Check pmd.adj 
	if( all(is.na(pmd.adj)) ) {
		pmd.adj.tmp <- partial.mahalanobis(object@x, object@mu, object@S)
		pmd.adj <- pmd.adj.tmp@pmd.adj
	}
	return(pmd.adj)
})
setGeneric("getDist", function(object) standardGeneric("getDist"))
setMethod("getDist", "CovRobMiss", function(object){
	pmd <- object@pmd
	## Check pmd.adj 
	if( all(is.na(pmd)) ) {
		pmd.adj.tmp <- partial.mahalanobis(object@x, object@mu, object@S)
		pmd <- pmd.adj.tmp@pmd
	}
	return(pmd)
})
setGeneric("getDim", function(object) standardGeneric("getDim"))
setMethod("getDim", "CovRobMiss", function(object) object@pu)
setGeneric("getLocation", function(object) standardGeneric("getLocation"))
setMethod("getLocation", "CovRobMiss", function(object) object@mu)
setGeneric("getScatter", function(object) standardGeneric("getScatter"))
setMethod("getScatter", "CovRobMiss", function(object) object@S)
setGeneric("getMissing", function(object) standardGeneric("getMissing"))
setMethod("getMissing", "CovRobMiss", function(object)which( object@pu == 0))
setGeneric("getOutliers", function(object, cutoff) standardGeneric("getOutliers"))
setMethod("getOutliers", "CovRobMiss", function(object, cutoff){
	pmd.adj <- object@pmd.adj
	## Check pmd.adj 
	if( all(is.na(pmd.adj)) ) {
		pmd.adj.tmp <- partial.mahalanobis(object@x, object@mu, object@S)
		pmd.adj <- pmd.adj.tmp@pmd.adj
	}
	if(missing(cutoff)) cutoff <- 0.99
	threshold <- qchisq( cutoff, df=object@p)
	return( which( pmd.adj > threshold) )
})

## S4 objects specific to CovRobMissSc objects extends to GSE and emve
setGeneric("getScale", function(object) standardGeneric("getScale"))
setMethod("getScale", "CovRobMissSc", function(object) object@sc)



## S4 method of plotting
setGeneric("plot")
setMethod("plot", signature(x="CovRobMiss", y="missing"), function(x, y="missing",
                                which=c("all", "index", "qqchisq", "dd"),
                                ask = (which=="all" && dev.interactive(TRUE)),
                                cutoff = 0.99, 
								xlog10 = FALSE, ylog10 = FALSE)
{
	pmd.adj <- x@pmd.adj
	## Check pmd.adj 
	if( all(is.na(pmd.adj)) ) {
		pmd.adj.tmp <- partial.mahalanobis(x@x, x@mu, x@S)
		pmd.adj <- pmd.adj.tmp@pmd.adj
	}
	obj <- list(pmd.adj=pmd.adj, p=x@p, x=x@x)

    which <- match.arg(which)

    op <- if (ask) par(ask = TRUE) else list()
    on.exit(par(op))
	
    ## qq-plot of the mahalanobis distances versus the
    ## quantiles of the chi-squared distribution
    if(which == "all" || which == "qqchisq") {
        print(.qqplot.pmdadj( obj, cutoff, xlog10, ylog10 ) )
    }

    ## index plot of partial square mahalanobis distances
    if(which == "all" || which == "index") {
        print(.distplot.pmdadj( obj, cutoff, ylog10 ))
    }
	
	## dd plot of partial square mahalanobis distances
    if(which == "all" || which == "dd") {
        print(.ddplot.pmdadj( obj, cutoff, xlog10, ylog10 ))
    }
}) 




