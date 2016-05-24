
# 
# Ingmar Visser, 11-6-2008
# 

# Changes
# - added lin.upper and lin.lower slots to these objects

# 
# MIX.FITTED CLASS
# 

setClass("mix.fitted",
	representation(message="character", # convergence information
		conMat="matrix", # constraint matrix on the parameters for general linear constraints
		lin.upper="numeric", # upper bounds for linear constraint
		lin.lower="numeric", # lower bounds for linear constraints
		posterior="data.frame" # posterior probabilities for the states
	),
	contains="mix"
)

setClass("mix.fitted.classLik",
	representation(
	    message="character", # convergence information
		conMat="matrix", # constraint matrix on the parameters for general linear constraints
		lin.upper="numeric", # upper bounds for linear constraint
		lin.lower="numeric", # lower bounds for linear constraints
		posterior="data.frame" # posterior probabilities for the states
	),
	contains="mix"
)

# accessor functions

setMethod("posterior","mix.fitted",
	function(object) {
		return(object@posterior)
	}
)

setMethod("show","mix.fitted",
	function(object) {
		cat("Convergence info:",object@message,"\n")
		print(logLik(object))
		cat("AIC: ", AIC(object),"\n")
		cat("BIC: ", BIC(object),"\n")
	}
)


# 
# Ingmar Visser, 23-3-2008
# 

# 
# DEPMIX.FITTED CLASS
# 

setClass("depmix.fitted",
	representation(message="character", # convergence information
		conMat="matrix", # constraint matrix on the parameters for general linear constraints
		lin.upper="numeric", # upper bounds for linear constraints
		lin.lower="numeric", # lower bounds for linear constraints
		posterior="data.frame" # posterior probabilities for the states
	),
	contains="depmix"
)

setClass("depmix.fitted.classLik",
	representation(message="character", # convergence information
		conMat="matrix", # constraint matrix on the parameters for general linear constraints
		lin.upper="numeric", # upper bounds for linear constraints
		lin.lower="numeric", # lower bounds for linear constraints
		posterior="data.frame" # posterior probabilities for the states
	),
	contains="depmix"
)
# accessor functions

setMethod("posterior","depmix.fitted",
	function(object) {
		return(object@posterior)
	}
)

setMethod("show","depmix.fitted",
	function(object) {
		cat("Convergence info:",object@message,"\n")
		print(logLik(object))
		cat("AIC: ", AIC(object),"\n")
		cat("BIC: ", BIC(object),"\n")
	}
)
