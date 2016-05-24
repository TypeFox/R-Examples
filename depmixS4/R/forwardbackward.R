# 
# Ingmar Visser
# 
# FORWARD-BACKWARD function, user interface, 10-06-2008
# 

setMethod("forwardbackward","depmix",
	function(object, return.all=TRUE, useC=TRUE, ...) {
		fb(init=object@init,A=object@trDens,B=object@dens,ntimes=ntimes(object), 
			homogeneous=object@homogeneous,return.all=return.all,useC=useC)
	}
)

setMethod("forwardbackward","mix",
	function(object, return.all=TRUE, useC=TRUE, ...) {
		fb(init=object@init,matrix(0,1,1),B=object@dens,ntimes=ntimes(object), 
			homogeneous=TRUE,return.all=return.all,useC=useC)
	}
)


