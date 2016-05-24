.validRSeed <- function(object) {
	return(TRUE)
}
setClass("RSeed",
			representation(
			env="environment"
			),
			validity=.validRSeed
)

setMethod("show", signature(object = "RSeed"),
	function(object){
		cat("==Original Model:\n")
		print(model_original(object))
		cat("==Connected Component Cut Off:\n")
		print(connectedComponentCutOff(object))
		cat("==currencyMetabolites:\n")
		print(currencyMetabolites(object))
		cat("\n==Edited Model:\n")
		print(model_edited(object))
		cat("\n==Model Changes:\n")
		if(!is.null(model_changes(object))){
			print(sapply(model_changes(object), length))
		}
		else{
			cat(" - \n")
		}
		cat("\n==Graph of the Network:\n")
		print(graph_network(object))
		
		cat("\n==List of the Strong Connected Components:\n")
		if(!is.null(scc(object))){
			print(length(scc(object)))
			print(summary(sapply(scc(object), length)))
		}
		else{
			cat(" - \n")
		}
		
		cat("\n==Graph of the Strong Connected Components:\n")
		print(graph_scc(object))
		cat("\n==List of the Source Compounds:\n")
		
		if(!is.null(scc_sizes(object))){
			print(length(list_sc(object)))
			print(summary(scc_sizes(object)[list_sc(object)]))
		}
		else{
			cat(" - \n")
		}
		
	}
)



slots <- c("model_original", "model_edited", "model_changes", "graph_network", "graph_scc", "scc", "scc_sizes", "list_sc", "connectedComponentCutOff", "currencyMetabolites")

for( i in slots){

	eval(parse(text=paste("\
		setGeneric(i,\
				function (rs) standardGeneric(\"", i, "\")\
		)\
		setGeneric(\"", i, "<-\",\
				function (rs, value) standardGeneric(\"", i, "<-", "\")\
		)\
		", 
		
		"setMethod(\"", i, "\", signature(rs = \"RSeed\"),\
			function(rs){\
				if(!exists(envir=rs@env, x=\"", i, "\", inherits=FALSE)){\
					return(NULL)\
				}\
				else{\
					return(get(\"", i, "\", rs@env, inherits=FALSE))\
				}\
			}\
		)\
		setReplaceMethod(\"", i, "\", signature(rs = \"RSeed\"),\
			function(rs, value){\
				assign(\"", i, "\", value, envir=rs@env, inherits=FALSE)\
				return(rs)\
			}\
		)",
		sep="", collapse=""))
	)
}

setGeneric("model_used",
		function (rs) standardGeneric("model_used")
)

setMethod("model_used", signature(rs = "RSeed"),
	function(rs){
		if(!is.null(model_edited(rs))){
			return(model_edited(rs))
		}
		else{
			return(model_original(rs))
		}
	}
)


rm(slots)































