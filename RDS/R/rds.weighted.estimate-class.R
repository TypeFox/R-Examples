################################################################################
# Class definition for "rds.weighted.estimate".   This class will be returned as the
# result of computing RDS estimates with no confidence interval.    
################################################################################

setClass("rds.weighted.estimate",
		representation(estimate="numeric",
				weights="numeric",
				outcome.variable="character",
				mean.group.weights="numeric",
				weight.type="character",
				subset="logical",
				csubset="character"))


setMethod("initialize",
		"rds.weighted.estimate",
		function(.Object,
				estimate,
				weights,
				outcome.variable,
				weight.type,
				subset,
				csubset) {
			
			.Object@estimate <- estimate
			.Object@weights <- weights
			.Object@outcome.variable <- outcome.variable
			.Object@weight.type <- weight.type
			se <- substitute(subset)
			if(is.null(se))
				subset <- rep(TRUE,length(weights))
			.Object@subset <- subset
			.Object@csubset <- csubset
			return(.Object)
		})

# A "friendly" constructor.  

rds.weighted.estimate <- function(estimate, weights,
		outcome.variable,weight.type,subset,csubset){
	new("rds.weighted.estimate", estimate, weights,
			outcome.variable,weight.type,subset,csubset)
}


setMethod("print",signature="rds.weighted.estimate",
		definition=function(x){
			cat(x@weight.type,"Estimate for",x@outcome.variable,
	switch(((x@csubset=="")|(x@csubset=="NULL"))+1,paste("[",x@csubset,"]",sep=""),NULL),
			"\n")
			if(!is.null(attr(x@estimate,"EL.se"))){
				se <- attr(x@estimate,"EL.se")
				attr(x@estimate,"EL.se") <- NULL
				print(x@estimate)
				if(se > 0){
					cat("Estimated standard error","\n")
					print(se)
				}
			}else{ 
				print(x@estimate)
			}
		}
)

setMethod("show",signature="rds.weighted.estimate",
		definition=function(object){
			cat(object@weight.type,"Estimate for",object@outcome.variable,
	switch(((object@csubset=="")|(object@csubset=="NULL"))+1,paste("[",object@csubset,"]",sep=""),NULL),
				"\n")
			if(!is.null(attr(object@estimate,"EL.se"))){
				se <- attr(object@estimate,"EL.se")
				attr(object@estimate,"EL.se") <- NULL
				print(object@estimate)
				if(se > 0){
					cat("Estimated standard error","\n")
					print(se)
				}
			}else{ 
				print(object@estimate)
			}
		}
)




