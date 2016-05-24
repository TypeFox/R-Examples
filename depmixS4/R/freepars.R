# depends on getpars(object) & getdf
setMethod("freepars","mix",
	function(object) {
		free <- getdf(object@prior)
		ns <- nstates(object)
		nresp <- nresp(object)
		for(i in 1:ns) {
			for(j in 1: nresp) {
				free <- free + getdf(object@response[[i]][[j]])
			}
		}		
		free
	}
)

# depends on getpars(object) & getdf
setMethod("freepars","depmix",
	function(object) {
		free <- getdf(object@prior)
		ns <- nstates(object)
		nresp <- nresp(object)
		for(i in 1:ns) {
			free <- free + getdf(object@transition[[i]])
			for(j in 1: nresp) {
				free <- free + getdf(object@response[[i]][[j]])
			}
		}
		free
	}
)

# depends on nlin(object) and getpars(object)
setMethod("freepars","depmix.fitted",
	function(object) {
		free <- sum(!getpars(object,which="fixed"))
 		free <- free-nlin(object) 
		free
	}
)

# depends on nlin(object) and getpars(object)
setMethod("freepars","mix.fitted",
	function(object) {
		free <- sum(!getpars(object,which="fixed"))
		free <- free-nlin(object) 
		free
	}
)

setMethod("nlin","mix.fitted",
	function(object) {
		conMat <- object@conMat[which(object@lin.lower==object@lin.upper),,drop=FALSE]
		if(nrow(conMat)==0) nlin <- 0 
		else nlin <- qr(conMat)$rank
		nlin
	}
)

setMethod("nlin","depmix.fitted",
	function(object) {
		conMat <- object@conMat[which(object@lin.lower==object@lin.upper),,drop=FALSE]
		if(nrow(conMat)==0) nlin <- 0 
		else nlin <- qr(conMat)$rank
		nlin
	}
)

