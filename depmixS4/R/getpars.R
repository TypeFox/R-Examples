setMethod("getpars","mix",
	function(object,which="pars",...) {
		parameters <- getpars(object@prior,which=which)
		for(i in 1:object@nstates) {
			for(j in 1:object@nresp) {
				parameters <- c(parameters,getpars(object@response[[i]][[j]],which=which))
			}
		}
		return(parameters)
	}
)

setMethod("getpars","depmix",
	function(object,which="pars",...) {
		parameters <- getpars(object@prior,which=which)
		for(i in 1:object@nstates) {
			parameters <- c(parameters,getpars(object@transition[[i]],which=which))
		}
		for(i in 1:object@nstates) {
			for(j in 1:object@nresp) {
				parameters <- c(parameters,getpars(object@response[[i]][[j]],which=which))
			}
		}
		return(parameters)
	}
)