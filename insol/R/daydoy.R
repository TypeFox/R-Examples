if ( !isGeneric("daydoy") ) {
	setGeneric("daydoy", function(x, month, day)
		standardGeneric("daydoy"))
}

setMethod('daydoy', signature(x="missing" ), 
	function(x){
		cat("USAGE: daydoy(year,month,day) \n")
        return()
	}
)

setMethod('daydoy', signature(x="numeric" ),
	function(x,month,day){ 
		return(as.POSIXlt(ISOdate(x,month,day))$yday+1)
	}
)

setMethod('daydoy', signature(x="POSIXct" ), 
	function(x){
		return(as.POSIXlt(x)$yday+1)
	}
)