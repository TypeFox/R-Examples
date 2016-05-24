setMethod("hiplarSet", signature(var = "character", val = "numeric"),
	function(var, val) .Call("hiplarSet", var, val, PACKAGE="HiPLARM")
)

hiplarShow <- function() .Call("hiplarShow", PACKAGE="HiPLARM")
