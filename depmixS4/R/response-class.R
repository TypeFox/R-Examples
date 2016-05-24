
# 
# response, transition and prior models for DEPMIX models
# 23-3-2008
# Maarten Speekenbrink & Ingmar Visser
# 

# 
# RESPONSE CLASS
# 

setClass("response",
	representation(parameters="list",
		fixed="logical",
		y = "matrix",
		x = "matrix",
		npar = "numeric", # this is not really needed as it simply is length(unlist(parameters))
		constr = "ANY"
	)
)

#
# RESPONSE CLASS METHODS
# 

setMethod("npar","response",
	function(object) {
		return(object@npar)
	}
)

setMethod("getdf","response",
	function(object) {
		return(sum(!object@fixed))
	}
)

