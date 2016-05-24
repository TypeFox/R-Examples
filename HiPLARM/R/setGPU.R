### Function to set our GPU Flag

setMethod("setGPU", signature(proc = "logical", interface = "character"),
	function(proc, interface) .Call("setGPU", proc, interface, PACKAGE="HiPLARM")
)


setMethod("setGPU", signature(proc = "numeric", interface = "character"),
	function(proc, interface) .Call("setGPU", proc, interface, PACKAGE="HiPLARM")
)

setMethod("setGPU", signature(proc = "numeric", interface = "missing"),
	function(proc, interface= "CPU") .Call("setGPU", proc, interface, PACKAGE="HiPLARM")
)

setMethod("setGPU", signature(proc = "logical", interface = "missing"),
	function(proc, interface = "CPU") .Call("setGPU", proc, interface = "CPU", PACKAGE="HiPLARM")
)

#setMethod("setinterface", signature(proc = "character"),
#	function(proc) .Call("interfaceSet", proc, PACKAGE="HiPLARM")
#)

#setGPU <- function(proc){ .C("GPUset", as.integer(proc), PACKAGE="HiPLARM") }


